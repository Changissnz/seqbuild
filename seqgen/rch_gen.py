"""
file contains code to extend the <RChainHead> class from the 
library <morebs2>
"""
from seqbuild.poly_output_fitter_ import * 
from morebs2.matrix_methods import is_vector
from morebs2.relevance_functions import RCInst,RChainHead

MRIF_VARMAP = {CEPoly:"v",\
    LinCombo:"x"}

def partitioned_vecmul(v1,v2): 
    return -1

def partitioned_vecdot(v1,v2):
    return -1

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
            self.acc_queue = acc_queue 
            self.qcap = queue_capacity

            self.mutgen = [{} for _ in range(len(self.rch.s))] 
            self.update_log = {}
            self.ctr = 0
            return 

    """
    main method
    """
    def apply(self,x): 
        assert is_vector(x) or type(x) in \
            {int,np.int32,np.int64}

        self.rch.apply(x)
        vx = deepcopy(self.rch.vpath)

        for v in self.rch.vpath: 
            if type(v) != np.ndarray: 
                self.acc_queue.append(v)
            else: 
                self.acc_queue.extend(v) 

        self.ctr += 1 
        self.update()
        return 

    #-------------------------------------

    def auto_add_mutable(self):
        return -1 

    def add_mutable(self,mg,rci_index,var_idn): 
        assert type(mg) == MutableRInstFunction
        assert var_idn in {'cf','dm','r'} 
        assert rci_index < len(self.mutgen) and \
            rci_index >= 0
        self.mutgen[rci_index][var_idn] = mg

    #--------------------- methods for updating 

    def update(self): 

        ml = self.mutable2update_list()

        for ml_ in ml: 
            self.update_idn(ml_[0],ml_[1])
        return 

    def update_idn(self,rci_index,var_idn):
        assert var_idn in {'r','cf','dm'}

        q = self.fetch_varlist_for_idn(rci_index,var_idn)

        # case: update reference value
        if var_idn == 'r': 
            self.rch.load_update_vars(q)
            self.tmpset_rch_updatepath(rci_index,len(q))
            self.rch[rci_index].inst_update_var() 

        # case: update function 
        else: 
            self.mut_gen[rci_index][var_idn].update(q)
            fx = self.mut_gen[rci_index][var_idn].apply
            self.rch[rci_index].update_var(var_idn,fx)
        return 

    def tmpset_rch_updatepath(self,rci_index,n): 
        x = [i for i in range(n)]
        self.rch.updatePath = {rci_index: x}
        return

    def fetch_varlist_for_idn(self,rci_index,var_idn):
        if len(self.acc_queue) == 0: 
            print("[!] none in queue for op.")
            return 

        vl = self.mutgen[rci_index][var_idn]
        d = vl.dim()
        return self.choose_n(d) 

    def choose_n(self,n):
        l = len(self.acc_queue)
        assert l > 0 
        q = []
        while n > 0:  
            i = self.prg() % l
            q.append(self.acc_queue[i])
            n -= 1
        return q 

    def mutable2update_list(self): 
        q = [] 
        for (i,x) in enumerate(self.mutgen): 
            for k,v in x.items():
                y = self.ctr % v.update_freq
                if y == 0:
                    q.append((i,k))
        return q 

    def __next__(self):
        return -1

    def change_dim(self): 
        return -1 