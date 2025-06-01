"""
file contains code to extend the <RChainHead> class from the 
library <morebs2>
"""
from seqbuild.poly_output_fitter_ import * 

MRIF_VARMAP = {CEPoly:"v",\
    LinCombo:"x"}

class MutableRInstFunction:

    def __init__(self,base_func,update_freq:int): 
        assert type(base_func) in MRIF_VARMAP
        assert type(update_freq) == int and update_freq >= 0 
        self.base_func = base_func 
        self.update_freq = update_freq 
        return

    def dim(self):
        q = getattr(self.base_func,MRIF_VARMAP[type(self.base_func)])
        return len(q) 

    def update(self,new_var): 
        assert type(new_var) == np.ndarray
        assert len(new_var) == self.dim()
        setattr(self.base_func,MRIF_VARMAP[type(self.base_func)],new_var)
        return

    def apply(self,x): 
        return self.base_func.apply(x) 

class RCHAccuGen: 

    def __init__(self,rch,perm_prg,acc_queue=[],\
        queue_capacity:int=1000):
            assert type(acc_queue) == list 
            assert type(queue_capacity) == int and queue_capacity > 1
            self.rch = rch 
            self.perm_prg = perm_prg 
            self.acc_queue = [] 
            self.update_log = {}
            return 

        def fetch_varlist(self):
            return -1 

        def apply(self,x): 
            self.rch.apply(x)
            vx = deepcopy(self.rch.vpath)
            self.acc_queue.extend(vx)
            self.update()
            return 

        def update(self): 
            return -1

        def __next__(self):
            return -1