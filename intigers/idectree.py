from .pof_autogen import * 
from .intfactor import * 

"""
Node structure used as unit for a decision tree constructed 
from an <IntSeq>. Structure is associated with functions, 
an input (operation) function `entryf` for an incoming value 
transmitted by a parent and an output (classifier) function `travf` 
to map the incoming value to the next node. 
"""
class IDecTNode: 

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
converts an <IntSeq> instance to a tree (directed graph) T. 
T satisfies a leaf requirement `l` XOR  a depth requirement `d`. 
"""
class IntSeqToTree: 

    def __init__(self,intseq,l:int,d:int): 
        assert type(intseq) == IntSeq
        assert type(l) == type(None) or type(d) == type(None)
        self.intseq = IntSeq 
        self.l = l 
        self.d = d
        self.leaf_first = type(l) != type(None) 
        if self.leaf_first: 
            assert type(self.l) == int and self.l > 0
            assert self.d <= len(self.intseq) - 1 
        else: 
            assert type(self.d) == int and self.d >= 0
            assert self.d <= len(self.intseq) - 1 

    def init_root(self): 
        return -1 

    def init_lf(self): 
        return -1 

    def init_df(self):
        return -1 

    #---------------------- splitter functions 

    def split_node(self,node): 
        return -1 

    def factor_split(self,S,cat_size:int): 
        assert len(S) > 0
        assert cat_size >= 1 
        s = ISFactorSetOps(list(S))

        cat2size = {}
        s.factor_count_()
        return -1 

    def poly_split(self,S,poly_size:int): 
        return -1 