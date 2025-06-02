"""
file contains code to extend the <RChainHead> class from the 
library <morebs2>
"""
from seqbuild.poly_output_fitter_ import * 
from morebs2.matrix_methods import is_vector
from morebs2.relevance_functions import RCInst,RChainHead
from morebs2.measures import zero_div 

MRIF_VARMAP = {CEPoly:"v",\
    LinCombo:"x"}

DM_FUNC_LIST = [np.dot,mul,safe_div,add,sub]

zero_div0 = lambda num,denum: zero_div(num,denum,0)

def safe_div(V1,V2):
    stat1 = is_vector(V1)
    stat2 = is_vector(V2) 

    if not stat1 and not stat2: 
        return zero_div0(V1,V2)

    if stat1 and stat2:
        assert len(V1) == len(V2)
        q = []
        for x,x2 in zip(V1,V2):
            q2 = zero_div0(x,x2)
            q.append(q2)
        return np.array(q)

    i = 0
    VX = None
    f = None 

    if stat1:
        VX = V1
        f = lambda x: zero_div0(x,V2)
    else: 
        VX = V2
        f = lambda x: zero_div0(V1,x) 

    q = [] 
    for x in VX:
        q.append(f(x))
    return q 
    
def safe_npint32_value(v):
    r = (np.iinfo(np.int32).min,np.iinfo(np.int32).max)
    v_ = modulo_in_range(v,r) 
    return np.int32(v_)

def safe_npint32_vec(V):
    return np.array([safe_npint32_value(v_) for \
        v_ in V],dtype=np.int32) 

def safe_npint32__prg_vec(prg,sz):
    v = np.zeros((sz,),dtype=np.int32) 
    for i in range(sz):
        v[i] = safe_npint32_value(prg()) 
    return v 

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
            self.update_log = defaultdict(defaultdict)
            self.ctr = 0
            return 

    @staticmethod
    def one_new_instance(prg,mutrate=0.5): 
        return -1 

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
        assert var_idn in {'cf','rf'} 
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
        assert var_idn in {'rf','cf'}

        q = self.fetch_varlist_for_idn(rci_index,var_idn)

        # case: update reference value
        if var_idn == 'rf': 
            self.rch.load_update_vars(q)
            self.tmpset_rch_updatepath(rci_index,len(q))
            self.rch[rci_index].inst_update_var() 

        # case: update function 
        else: 
            self.mut_gen[rci_index][var_idn].update(q)
            fx = self.mut_gen[rci_index][var_idn].apply
            self.rch[rci_index].update_var(var_idn,fx)

        if rci_index not in self.update_log: 
            self.update_log[rci_index] = defaultdict(int) 
        self.update_log[rci_index][var_idn] += 1 
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