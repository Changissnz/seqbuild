from intigers.mdr_v2 import * 
from intigers.tvec import * 
from intigers.intfactor import * 
from intigers.extraneous import DEFAULT_PAIRWISE_OPS
from mini_dm.iseq import * 
from mini_dm.nsfr import * 
from intigers.extraneous import to_trinary_relation,to_trinary_relation_v2,zero_div0
from morebs2.matrix_methods import cr 
from collections import deque 
from morebs2.numerical_generator import sign_preserving_modulo,prg_decimal
from .ssi_load import OpTriGenLite

DEFAULT_SHADOW_FITTERS = {"mdr","mdrv2","tvec","fvec","optri"}
DEFAULT_SHADOWGEN_BASESEQ_LENGTH_RANGE = [7,91]
DEFAULT_SHADOWGEN_NUM_SHADOWS_RANGE = [2,8]
DEFAULT_SHADOWGEN_VECTOR_CACHE_SIZE = 52543
DEFAULT_SHADOWGEN_MDR_MAX_ABSMULT = 44 
DEFAULT_SHADOWGEN_FVEC_MAX_ABSMULT = 59

DEFAULT_SHADOWGEN_MAX_ABS = ceil(9.97 ** 6) 

"""
Shadown Pseudo-Random Number Generator. 

Named because of "shadowing" process over cyclical streaming of `file_path` containing a numerical 
sequence. There are these 5 types of fitters: `mdr`,`mdrv2`,`tvec`,`fvec`,`optri`. 

For every one of these five fitters, uses an auxiliary PRNG `prg` to make decisions that add noise 
to the underlying form of some (x in DEFAULT_SHADOWGEN_BASESEQ_LENGTH_RANGE) numbers from the 
`file_path` stream. With the parameters of the underlying form changed, this PRNG then outputs a 
typically differing number sequence from that of `base_vec`. 

For `mdr` and `mdrv2` (modulo decomposition representation), the underlying form is the segments of 
(multiple,additive,modulo) connecting the contiguous elements of the current `base_vec`. 

For `tvec`, the ternary vector of the current `base_vec` is changed. 

For `fvec`, the factors of every element in `base_vec` are changed into a new factor set F, and the 
resulting output element from that element is the cumulative product of F. 

For `optri`, this PRNG selects a forward function and backward function from (+,-,*,/). Then it 
forms a matrix of dimension |base_vec| x |base_vec|, roughly half of the elements comprised from the 
forward function and the remainder from the backward function. This matrix is flattened into a sequence 
and then shuffled using the `prg`. The resultant is the expressed output values from the `base_vec`. 

NOTE: `fvec`, due to the process collecting factor sets for every element of the current `base_vec`, 
    runs considerably slower than the other `fitting_struct`s, on a scale of at least 50X. 
"""
class ShadowGen:

    def __init__(self,prg,file_path,fitting_struct,cast_func = cr):
        assert type(prg) in {MethodType,FunctionType}
        assert os.path.exists(file_path)
        assert fitting_struct in DEFAULT_SHADOW_FITTERS

        self.prg = prg 
        self.file_path = file_path

        file_obj = open(self.file_path,'r') 
        is_periodic = True 
        self.nsfr = NSFileReader(file_obj,cast_func,is_periodic)
        self.fitting_struct = fitting_struct 
        self.cast_func = cast_func 

        self.fitter = None 

        self.queue = deque() 
        self.vec_cache = deque()  
        self.pv_size = 0 

        self.base_vec = None 
        return

    def __next__(self): 

        if len(self.queue) == 0: 
            self.load_next_fitter() 

        r = self.queue.popleft() 
        return self.cast_func(\
            sign_preserving_modulo(r,DEFAULT_SHADOWGEN_MAX_ABS)) 

    def load_next_fitter(self): 

        self.load_next_base_vec()

        V = None 
        if self.fitting_struct in {"mdr","mdrv2"}:  
            V = self.mdr_derivation() 
        elif self.fitting_struct == "tvec": 
            V = self.tvec_derivation() 
        elif self.fitting_struct == "fvec": 
            V = self.fvec_derivation() 
        elif self.fitting_struct == "optri":
            V = self.optri_derivation() 
        else: 
            assert False 

        self.queue.extend(V) 
        self.update_vec_cache(V) 

    def update_vec_cache(self,V): 
        q = len(V) 
        self.pv_size += q 

        V = np.array([sign_preserving_modulo(v,DEFAULT_SHADOWGEN_MAX_ABS) for v in V])
        self.vec_cache.append(V) 

        while self.pv_size > DEFAULT_SHADOWGEN_VECTOR_CACHE_SIZE: 
            x = self.vec_cache.popleft() 
            self.pv_size -= len(x)

    def close(self): 
        self.nsfr.close() 
        del self 

    #---------------------------------------- methods for loading up next base vector 

    def load_next_base_vec(self): 
        l = modulo_in_range(int(self.prg()),DEFAULT_SHADOWGEN_BASESEQ_LENGTH_RANGE)
        self.base_vec = np.array([next(self.nsfr) for _ in range(l)])
        self.add_shadows_to_base_vec() 

    def add_shadows_to_base_vec(self): 
        l_ = len(self.vec_cache)
        if l_ == 0: 
            return 

        l = modulo_in_range(int(self.prg()),DEFAULT_SHADOWGEN_NUM_SHADOWS_RANGE) 
        l = min([l_,l])

        indices = [i for i in range(l_)] 
        prg_ = prg__single_to_int(self.prg)

        indices = prg_choose_n(indices,l,prg_,is_unique_picker=True)

        Q = [self.vec_cache[i] for i in indices]
        
        original_length = len(self.base_vec)
        for q in Q: 
            self.base_vec = modulated_vec_op(self.base_vec,q,add)
            self.base_vec = self.base_vec[:original_length] 
        return 

    #----------------------------------------------------- derivations on fitters for current base vector 

    def load_mdr(self): 

        V = [sign_preserving_modulo(v,DEFAULT_SHADOWGEN_MAX_ABS) for v in self.base_vec]  
        t = 1 
        if self.fitting_struct == "mdrv2":  
            md = ModuloDecompV2(IntSeq(V),False,DEFAULT_SHADOWGEN_MDR_MAX_ABSMULT) 
            t = 2 
        else: 
            md = ModuloDecomp(IntSeq(V),DEFAULT_SHADOWGEN_MDR_MAX_ABSMULT) 
            md.merge(False)  
        self.fitter = ModuloDecompRepr(md,t)
        return 

    def mdr_derivation(self): 
        self.load_mdr() 

        d = prg_decimal(self.prg,[0.,1.])

        # case: shift the partition 
        if d > 0.5: 
            lx = len(self.fitter.afs_prt) 
            lx = lx ** 2 
            s = int(round(self.prg() % lx))
            q = self.fitter.shift_afs_prt_(s) 
            self.fitter.afs_prt = q 
        # case: add noise to partition 
        else: 
            self.fitter.noise_to_afs_prt(self.prg,True)
        return self.fitter.reconstruct() 

    def tvec_derivation(self): 
        # get the current trinary vector 
        T = gleqvec(self.base_vec)

        # get a trinary vector different from this one 
        Q = list(generate_m_unique_trinary_vectors(len(T),2,self.prg,attempt_ratio=3.0)) 
        if np.all(T == Q[0]):
            T = Q[1] 
        else: 
            T = Q[0] 

        return ternary_adjustment(self.base_vec,T,self.prg)

    def fvec_derivation(self): 

        V = np.array(self.base_vec,dtype=int)
        isfso = ISFactorSetOps(V,int_limit=DEFAULT_INT_MAX_THRESHOLD,str_mode_full=True) 
        isfso.factor_count_() 
        factors = isfso.factors 

        L = []
        for (i,b) in enumerate(V): 
            # get factors for value 
            if len(factors[i]) == 0: 
                q = np.array([10,7,3]) 
            else: 
                q = np.array(sorted(factors[i])) 

            # add noise to factors by a unique vector 
            U = prg_unique_sequence(self.prg,len(q))
            q = q + np.array(U) 

            # cumulative product of q 
            L.append(np.cumprod(q)[-1]) 
        return np.array(L) 

    def optri_derivation(self): 
        L = IntSeq(self.base_vec)

        i = int(self.prg()) % len(DEFAULT_PAIRWISE_OPS) 
        j = int(self.prg()) % len(DEFAULT_PAIRWISE_OPS) 

        op1 = DEFAULT_PAIRWISE_OPS[i]
        op2 = DEFAULT_PAIRWISE_OPS[j]

        # make the square matrix of |base_vec| x |base_vec| 
        otl = OpTriGenLite(L,op1,op2) 
        V = list(otl.m.flatten()) 

        # shuffle the flattened matrix 
        V = prg_seqsort(V,self.prg) 
        return np.array(V) 