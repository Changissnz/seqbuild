from intigers.n2mvf_autogen import *
from morebs2.numerical_generator import prg__LCG,prg__single_to_nvec
from collections import deque 

DEFAULT_N2MV_GEN_FILLER_PRG2_MODULO = 7891

"""
The PRNG `prg` is the primary generator: this generator uses `prg` 
to produce input vectors of dimension n to output m-vectors from the 
current <N2MVectorFunction>. 

The main `__next__` function outputs values from its queue. 

When `queue` is empty, generator has its <N2MVectorFunctionGen> produce 
a <N2MVectorFunction> F. Then F outputs an m-vector. That m-vector is added 
to `queue`. After some k >= (s in `update_frequency_range`) output values, 
this generator instantiates another <N2MVectorFunctionGen>. 

See file<intigers.n2mvf_autogen> for more details on 
"""
class N2MVGen: 

    def __init__(self,nm_range,prg,prg2,update_frequency_range):
        assert is_valid_range(nm_range,True,False) 
        assert type(prg) in {FunctionType,MethodType} 
        assert type(prg2) in {FunctionType,MethodType,type(None)} 
        assert is_valid_range(update_frequency_range,True,False)

        self.nm_range = nm_range 
        self.prg = prg__single_to_int(prg) 

        self.prg2 = prg2 
        if type(self.prg2) != type(None):  
            self.prg2 = prg__single_to_int(prg2)   

        self.nvfg = None 

        self.queue = deque() 

        self.nvfg = None 
        self.current_input_dim = None 

        self.uf_range = update_frequency_range
        self.update_threshold = None 
        self.update_ctr = 0 
        self.filler_prg2() 
        self.next_NVMVFGen() 

    def __next__(self): 
        if len(self.queue) == 0: 
            self.next_vector() 
        return self.queue.popleft()

    def next_NVMVFGen(self): 

        n = modulo_in_range(int(self.prg()),self.nm_range) 
        m = modulo_in_range(int(self.prg()),self.nm_range) 

        nm = (n,m)
        self.current_input_dim = n 

        mode = "replace" if int(self.prg()) % 2 else "accumulate"
        self.nvfg = N2MVectorFunctionGen(nm,self.prg,self.prg2,mode)

        self.update_threshold = modulo_in_range(int(self.prg()),self.uf_range) 
        self.update_ctr = 0 

    def next_vector(self): 
        q = prg__single_to_nvec(self.prg,self.current_input_dim)() 
        F = self.nvfg.one_N2MVF() 
        mvec = F.apply(q)

        self.update_ctr += len(mvec)  

        self.queue.extend(mvec)  

        if self.update_ctr >= self.update_threshold: 
            self.next_NVMVFGen()

    """
    an LCG, in cases where `prg2` is null.
    """
    def filler_prg2(self): 
        if type(self.prg2) != type(None): 
            return 

        x = prg_unique_sequence(self.prg,4) 
        if x[-1] == 0: 
            x[-1] = DEFAULT_N2MV_GEN_FILLER_PRG2_MODULO
        self.prg2 = prg__LCG(x[0],x[1],x[2],x[3]) 

