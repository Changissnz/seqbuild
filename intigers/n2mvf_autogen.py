"""
file contains code that auto-generates a <N2MVectorFunction> instance. 
"""

from mini_dm.n2m_function import * 
from intigers.poly_output_fitter_ import * 
from intigers.mod_prng import prg__iterable
from intigers.extraneous import subvec
from morebs2.numerical_generator import prg__constant

"""
function input dimension; used specifically for 
<LinCombo> and <CEPoly> instances. 
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
shifts of the input n-vector. The dictionary `fmap`
is the mapping of abstract functions to indices, 
so class is not exclusively programmed for linear 
combinations and polynomial functions. 

The variables `index_shifter` and `subvec_size_shifter`
are parameter-less functions that return integers. They 
effectively determine each subvector of the input n-vector 
`x` that corresponds to each index in the output m-vector.
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

        # change the fmap from container of <CEPoly> + <LinCombo>
        # to that of abstract functions. 
        fmap_ = [] 
        for f in fmap: 
            if type(f) == CEPoly:
                fx = lambda x: f.apply(x[0]) 
            else: 
                fx = f.apply 
            fmap_.append(fx)

        # make the index shifter 
        index_shifter = prg__constant(1)

        return LCPVectorMap__TypeCShift(nm,fmap_,index_shifter,subvec_size_shifter)


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
        self.prg = prg2
        self.mode = mode 
