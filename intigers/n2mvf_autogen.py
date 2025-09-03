"""
file contains code that auto-generates a <N2MVectorFunction> instance.
There are no constraints on functions from n to m space, as long as 
said function `f` is capable of 
    `f(V_n) = V_m`, V_n and V_m are two subvectors associated with the 
            n-vector and m-vector, respectively. The specific sizes of 
            `V_n` and `V_m` are determined by user or machine input 
            for the corresponding (n,m) index-set pair. 

However, two popular mathematical forms are linear combinations 
and polynomials. These two forms are represented by 
    class<LCPVectorMap__TypeCShift>. 
Additionally, modulated vector operations are possible through 
class<ModulatedN2MVectorMap>. 
"""

from mini_dm.n2m_function import * 
from intigers.poly_output_fitter_ import * 
from intigers.mod_prng import prg__iterable
from mini_dm.ngram import *
from morebs2.numerical_generator import prg__constant
from intigers.prng_pw_op import * 

"""
outputs function input dimension; used specifically 
for <LinCombo> and <CEPoly> instances. 
"""
def lcp_dim(F): 
    assert type(F) in {LinCombo,CEPoly}

    if type(F) == CEPoly:
        return 1
    return F.degree()  

"""
generates functions that are either linear combinations 
or polynomial functions. Values generated use `prg` for 
input vector into <LinCombo> or <CEPoly> instances. The 
generator `bprg` is used to randomly select between 
<LinCombo> and <CEPoly>. 

These functions are generated in a sequence. Class is 
used for generating <LCPVectorMap__TypeCShift> instances. 
"""
class LCPGen: 

    def __init__(self,prg,bprg,subvec_size_shifter):
        assert type(prg) in {MethodType,FunctionType}
        assert type(bprg) in {MethodType,FunctionType}

        self.prg = prg 
        self.bprg = bprg
        self.subvec_size_shifter = subvec_size_shifter
        return 

    """
    main method
    """
    def generate_lcp_sequence(self,seq_size):
        fseq = []
        for i in range(seq_size):
            sz = self.subvec_size_shifter()
            assert sz >= 0 and type(sz) in {int,np.int32,np.int64}

            v = np.array([self.prg() for _ in range(sz)])
            q = int((self.bprg()) % 2) 
            fx = None 
            if q == 0: 
                fx = LinCombo(v)
            else: 
                v2 = np.arange(sz)[::-1] 
                v = np.array([v,v2]).T 
                fx = CEPoly(v)
            fseq.append(fx)
        return fseq

"""
linear-combination-and-polynomial vector map 
from n to m space. Map is capable of contiguous 
shifts of the input n-vector. The list of 
functions `fmap` is a 1-to-1 mapping, such that 
every j'th function corresponds to the m-vector's 
j'th index. Every function f is an abstract function 
operated by the call `f()`, so class is not exclusively 
programmed for linear combinations and polynomial functions. 
The variables `index_shifter` and `subvec_size_shifter`
are parameter-less functions that return integers. They 
effectively determine each contiguous subvector of the 
input n-vector `x` that corresponds to each index in the 
output m-vector.
"""
class LCPVectorMap__TypeCShift: 

    def __init__(self,nm,fmap,index_shifter,subvec_size_shifter):
        assert_nm(nm)
        assert type(fmap) == list 
        for x in fmap: assert type(x) in {MethodType,FunctionType}
        assert type(index_shifter) in {MethodType,FunctionType}
        assert type(subvec_size_shifter) in {MethodType,FunctionType}

        self.nm = nm 
        self.fmap = fmap
        self.index = 0
        self.index_shifter = index_shifter
        self.subvec_size_shifter = subvec_size_shifter
        return 

    def apply(self,x):
        assert len(x) == self.nm[0] 
        assert is_vector(x) 

        x2 = np.zeros((self.nm[1],),dtype=float) 

        for i in range(self.nm[1]): 
            sz = self.subvec_size_shifter()
            sv = np.array(subvec(x,self.index,sz))
            F = self.fmap[i]
            q = F(sv)
            x2[i] = q 
            self.index = (self.index + self.index_shifter()) % self.nm[0] 
        return x2

    # TODO: work on this.
    @staticmethod 
    def one_LCPVectorMap(nm,subvec_size_shifter,prg,prg2):
        lcpv_gen = LCPGen(prg,prg2,subvec_size_shifter)
        fmap = lcpv_gen.generate_lcp_sequence(nm[1]) 

        # get the dimension for the function sequence 
        dim_fmap = [lcp_dim(x) for x in fmap]
        subvec_size_shifter = prg__iterable(dim_fmap)

        def make_function(f0,is_poly:bool):
            if is_poly: 
                return lambda x: f0.apply(x[0]) 
            return lambda x: f0.apply(x) 

        # change the fmap from container of <CEPoly> + <LinCombo>
        # to that of abstract functions. 
        fmap_ = [] 
        for f in fmap: 
            fx = make_function(f,type(f) == CEPoly)
            fmap_.append(fx)

        # make the index shifter 
        index_shifter = prg__constant(1)

        return LCPVectorMap__TypeCShift(nm,fmap_,index_shifter,subvec_size_shifter)

"""
structure transforms the pairwise operator `modulated_vec_op` 
into a uni-wise operator. Uses `v` as the base vector (first 
argument to `modulated_vec_op`) and operator `op` for the 
application, 
        `f(x) = modulated_vec_op(v,x,op)`. 
Every output from `apply` is a vector of length `|v|`. 
"""
class ModulatedN2MVectorMap:

    def __init__(self,v,op):
        assert is_vector(v)
        assert len(v) > 0

        self.v = v
        self.op = op
        return
    
    def apply(self,x): 
        assert is_vector(x)
        q = modulated_vec_op(self.v,x,self.op)
        return q[:len(self.v)]

    @staticmethod 
    def one_instance(prg,l,op): 
        vx = [prg() for _ in range(l)]
        vx = np.array(vx) 
        return ModulatedN2MVectorMap(vx,op) 

DEFAULT_N2MVF_INDEX_DEGREE_RANGE = (2,20) 
DEFAULT_N2MVF_NSET_SIZE_RATIO_RANGE = (0.22,0.6)
DEFAULT_N2MVF_MSET_SIZE_RATIO_RANGE = (0.05,0.8)

class N2MVectorFunctionGen:

    def __init__(self,nm,prg,prg2,mode="replace"):
        assert_nm(nm) 
        assert type(prg) in {MethodType,FunctionType}
        assert type(prg2) in {MethodType,FunctionType}
        assert mode in {"replace","cumulative"}

        self.nm = nm  
        self.prg = prg 
        self.prg2 = prg2
        self.mode = mode 

    """
    main method 
    """
    def one_N2MVF(self,map_gen_attempts_per=100): 
        nset_min = int(ceil(DEFAULT_N2MVF_NSET_SIZE_RATIO_RANGE[0] * self.nm[0]))
        nset_max = int(ceil(DEFAULT_N2MVF_NSET_SIZE_RATIO_RANGE[1] * self.nm[0]))
        
        mset_min = int(ceil(DEFAULT_N2MVF_MSET_SIZE_RATIO_RANGE[0] * self.nm[1]))
        mset_max = int(ceil(DEFAULT_N2MVF_MSET_SIZE_RATIO_RANGE[1] * self.nm[1]))

        if nset_min == nset_max: nset_max = nset_max + 1
        if mset_min == mset_max: mset_max = mset_max + 1

        n2mgen = N2MIndexMapGen(self.nm,\
            DEFAULT_N2MVF_INDEX_DEGREE_RANGE,[nset_min,nset_max],\
            [mset_min,mset_max],self.prg,index_degree_is_geq=False)
        
        # max attempts is square that of `map_gen_attempts_per`
        n2mgen.make(num_iter=map_gen_attempts_per,\
            attempts_per_relation=map_gen_attempts_per) 
        M = list(n2mgen.map()) 
        M = sorted(M,key=lambda x: x[0] + x[1]) 
        M = prg_seqsort(M,self.prg) 

        fmap = [] 
        for m in M: 
            is_lcpvm = bool(int(self.prg2() % 2)) 

            n0 = len(string_to_vector(m[0])) 
            m1 = len(string_to_vector(m[1])) 

            fx = self.nm_to_vmap(n0,m1,is_lcpvm).apply 
            fmap.append(fx) 
            
        fmap = {fx: [i] for (i,fx) in enumerate(fmap)}
        return N2MVectorFunction(self.nm,M,fmap,self.mode)

    """
    generates a <LCPVectorMap__TypeCShift> or <ModulatedN2MVectorMap>
    instance, based on reference to `n` and `m`, the sizes for the 
    n-set and m-set that this method's output function corresponds to. 

    NOTE: algorithm always uses weighted pairwise operator for 
          <ModulatedN2MVectorMap>. 
    """
    def nm_to_vmap(self,n,m,is_lcpvm:bool):
        
        # case: LCPVectorMap__TypeCShift
        if is_lcpvm:
            subvec_size_shifter = self.one_subvec_size_shifter(n) 
            return LCPVectorMap__TypeCShift.one_LCPVectorMap(\
                (n,m),subvec_size_shifter,self.prg,self.prg2)

        # case: ModulatedN2MVectorMap
        op = prg__one_weighted_pairwise_operator(self.prg,\
            deepcopy(DEFAULT_PAIRWISE_OPS),deepcopy(DEFAULT_PAIRWISE_OPS),None)
        return ModulatedN2MVectorMap.one_instance(self.prg,m,op)

    """
    generates an iterable over a sequence of 
    subvector sizes, each no greater than 
    `n0`. 
    """
    def one_subvec_size_shifter(self,n0,sz=5): 
        s = [int(modulo_in_range(self.prg(),DEFAULT_N2MVF_INDEX_DEGREE_RANGE)) \
            for _ in range(sz)]
        s2 = [] 
        for s_ in s: 
            if s_ > n0:
                q = max([3,n0+1])
                s1 = modulo_in_range(s_,[2,q]) 
            else: 
                s1 = s_ 
            s2.append(s1) 
        return prg__iterable(s2)