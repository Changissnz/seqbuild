from .pof_autogen import * 
from .intfactor import * 
from morebs2.numerical_generator import prg_seqsort_ties

"""
Node structure used as unit for a decision tree constructed 
from an <IntSeq>. Structure is associated with functions, 
an input (operation) function `entryf` for an incoming value 
transmitted by a parent and an output (classifier) function `travf` 
to map the incoming value to the next node. 
"""
class IDecNode: 

    def __init__(self,idn:int,entryf,travf,root_distance:int): 
        assert type(idn) == int  
        self.idn = idn 
        self.entryf = entryf
        self.travf = travf 
        self.children = [] 
        self.rd = root_distance 
        self.reachability_key = None 
        self.con_key = set() 
        self.acc_queue = []
        return

    #----------------------- var-setter functions 

    def set_travf(self,travf):
        self.travf = travf 

    def set_rk(self,rk): 
        self.reachability_key = rk

    def set_conkey(self,s): 
        assert type(s) == set 
        self.con_key = s

    def clear_conkey(self): 
        self.con_key.clear()  

    #----------------------- children-access functions 

    def add_children(self,idn_seq,entryf_seq): 
        assert len(idn_seq) == len(entryf_seq)
        for (i,x) in enumerate(idn_seq): 
            q = IDecNode(x,entryf_seq[i],None,self.rd+1)
            self.children.append(q)

    def index_of(self,idn): 
        q = self.children
        for (i,q_) in enumerate(q): 
            if q_.idn == idn: return i 
        return -1

    def fetch_conn(self,idn): 
        i = self.index_of(idn) 
        if i == -1: return False         
        return self.children[i]

    #----------------------- input-output function operation 

    def operate(self,v): 
        return self.entryf(v) 

    def classify(self,v): 
        return self.travf(v) 

"""
represents a factor-based or polynomial based boolean classifier 
"""
class IDecTrFunc: 

    def __init__(self,conditional_list,conditional_type):  
        assert conditional_type in {"poly","factor"}
        if conditional_type == "poly": 
            for x in conditional_list: 
                assert type(x) in {PolyOutputFitterVar1,PolyOutputFitterVar2}
        else: 
            for x in conditional_list: 
                assert type(x) in {int,np.int32,np.int64} 

        self.cl = conditional_list
        self.ct = conditional_type
        return

    """
    boolean classification 
    """
    def bclassify(self,x):  
        f = self.factor_classify if self.ct == "factor" else \
            self.poly_classify

        for i in range(len(self.cl)): 
            if f(x,i): return 1 
        return 0

    def factor_classify(self,x,i): 
        assert self.ct == "factor" 
        return x / self.cl[i] == x // self.cl[i]

    def poly_classify(self,x,i): 
        assert self.ct == "poly" 
        if type(self.cl[i]) == PolyOutputFitterVar1:
            return self.cl[i].apply(x) == self.cl[i].c 
        return self.cl[i].apply(x) == self.cl[i].apply(self.cl[i].x1) 

"""
converts an <IntSeq> instance to a tree (directed graph) T. 
T satisfies a leaf requirement `l` XOR  a depth requirement `d`. 
"""
class IntSeq2Tree: 

    def __init__(self,intseq,l:int,d:int,prg): 
        assert type(intseq) == IntSeq
        assert type(l) == type(None) or type(d) == type(None)
        self.intseq = intseq 
        self.l = l 
        self.d = d
        self.prg = prg 
        self.leaf_first = type(l) != type(None) 
        if self.leaf_first: 
            assert type(self.l) == int and self.l > 0
            assert self.l <= len(self.intseq) - 1 
        else: 
            assert type(self.d) == int and self.d >= 0
            assert self.d <= len(self.intseq) - 1 

        self.factor_preproc()

    def factor_preproc(self):
        self.isfso = ISFactorSetOps(deepcopy(self.intseq.l),\
            int_limit=DEFAULT_INT_MAX_THRESHOLD)
        self.isfso.factor_count_()
        return 

    #----------------------- root initialization to satisfy depth or
    #----------------------- leaf requirement 

    def init_root(self): 
        return -1 

    def init_lf(self): 
        return -1 

    def init_df(self):
        return -1 

    def split_node(self,node): 
        return -1

    #---------------------- factor-splitter functions 

    """
    splits a sequence S of elements according to the `partition` (list) given. 

    Outputs a dictionary with key `partition index` and value as a sequence of  
    (factor, elements of S divisible by factor). 

    If `last_subset_isneg` is set to False, allocates an additional element for 
    the output. Otherwise, last partition is implicit (not part of output). 
    """
    def factor_split__partitioned(self,S,partition,last_subset_isneg:bool=True):
        assert type(partition) == list 
        assert len(partition) > 1 
        assert min(partition) > 0 
        assert sum(partition) == len(S) 

        self.fs = self.sort_factors_by_keys(S,orderng = "min",m=None) 
        # NOTE: ineff, method sorts a second time 
        self.fs = prg_seqsort_ties(self.fs,self.prg,vf=lambda x:x[1])

        px = partition if not last_subset_isneg else partition[:-1] 
        FS = dict() 
        for (i,p) in enumerate(px): 
            fx = self.factors_for_subset(S,p) 
            FS[i] = fx 
        return FS 

    """
    sequence of elements as 
    (factor,{elements of S divisible by factor}) 
    """
    def factors_for_subset(self,S,sz): 
        f = [] 
        S = set(S) 
        while sz > 0: 
            q = self.factor_for_size(sz,self.fs)
            if type(q) == type(None):
                raise ValueError("something wrong: num factors {} sz req {}".format(len(self.fs),sz))
                break 
            factor = self.fs.pop(q) 
            keys = self.isfso.factor_to_keys(factor[0]) 
            keys = keys.intersection(S) 
            S -= keys 
            self.update_sorted_factors(keys)
            self.isfso.remove_seq_elements(keys) 
            f.append((factor[0],keys))  
            sz -= len(keys)  
        return f 

    # TODO: ineff 
    def factor_for_size(self,sz,fs):
        if sz == 0:
            return None 

        if len(fs) == 0: 
            return None  

        index_range = [0,None]
        mindiff = float('inf') 
        for (i,fs_) in enumerate(fs):
            if fs_[1] > sz: 
                index_range[1] = i 
                break 

            if sz - fs_[1] < mindiff: 
                index_range[0] = i
                mindiff = sz - fs_[1]

        if type(index_range[1]) == type(None): 
            index_range[1] = len(fs) 

        q = modulo_in_range(self.prg(),(index_range[0],index_range[1]))
        return q 

    def sort_factors_by_keys(self,S,orderng = "min",m=None):
        assert orderng in {"min","max","medi"}
        assert len(S) > 0 

        if orderng in {"min","max"}: 
            index = 0 if orderng == "min" else -1
            r = self.isfso.dsort(pkeys=S)
            if orderng == "max": 
                r = r[::-1] 
        else: 
            median = 0.5 if type(m) == type(None) else m 
            r = self.isfso.median_sort(pkeys=S,r=m,fullpair_sequence=True)
        return r 

    def update_sorted_factors(self,elements_to_remove): 

        if len(self.fs) == 0: 
            return 

        for (i,x) in enumerate(self.fs): 
            c = 0 
            for etr in elements_to_remove:
                if etr / x[0] == etr // x[0]: 
                    c += 1
            x[1] -= c
        
        self.fs = [fs_ for fs_ in self.fs if fs_[1] > 0] 
        self.fs = sorted(self.fs,key=lambda x:x[1])

    # NOTE: alternative method to <factor_split__partitioned>
    # TODO: test 
    """
    chooses a factor by degree based on ordering of 'min','max','median'. 
    Used for conditional splits with integer sets. 
    """
    def csplitting_factor(self,S,orderng = "min",m=None):
        assert orderng in {"min","max","medi"}
        assert len(S) > 0 
        rx = self.sort_factors_by_keys(S,orderng,m) 
        return rx[0][0] 

    #-------------------------------------------------------------------

    def poly_split(self,S,poly_size:int): 
        return -1 
