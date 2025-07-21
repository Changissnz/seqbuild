"""
file contains code that auto-generates a <N2MVectorFunction> instance. 
"""

from mini_dm.n2m_function import * 
from intigers.poly_output_fitter_ import * 

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

    def generate_function_sequence(seq_size):
        fseq = []
        for i in range(seq_size):
            sz = self.subvec_size_shifter()
            v = np.array([self.prg() for _ in range(sz)])
            q = int((self.bprg) % 2) 
            fx = None 
            if q == 0: 
                fx = LinCombo(v).apply
            else: 
                v2 = np.arange(sz)[::-1] 
                v = np.array([v,v2]).T 
                fx = CEPoly(v).apply
            fseq.append(fx)
        return fseq

"""
linear-combination-and-polynomial vector map 
from n to m space. Map is capable of contiguous 
shifts of the input n-vector. The dictionary `fmap`
is the mapping of abstract functions to indices, 
so class is not exclusively programmed for linear 
combinations and polynomial functions. 
"""
class LCPVectorMap__TypeCShift: 

    def __init__(self,nm,fmap,index_shifter,subvec_size_shifter):
        assert_nm(nm)
        assert type(index_shifter) in {MethodType,FunctionType}
        assert type(subvec_size_shifter) in {MethodType,FunctionType}

        self.nm = nm 
        self.index = 0
        self.index_shifter = index_shifter
        self.subvec_size_shifter = subvec_size_shifter
        return 

    def fit(self,x):
        assert len(x) == self.nm[0] 
        assert is_vector(x) 

        x2 = np.zeros((self.nm[1],),dtype=float) 

        for i in range(self.nm[1]): 
            sz = self.subvec_size_shifter()
            sv = np.array(subvec(x,self.index,sz))
            F = self.index_to_function(i)
            q = F(sv)
            x2[i] = q 
            self.index = (self.index + self.index_shifter()) % self.nm[0] 
        return x2

    def index_to_function(self,i): 
        for k,v in self.fmap.items():
            if v == i:
                return k
        return None 

class ModulatedN2MVectorMap:

    def __init__(self):
        return -1 
        #modulated_vec_op(v1,v2,op)

DEFAULT_N2MVF_INDEX_DEGREE_RANGE = (1,20) 
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
