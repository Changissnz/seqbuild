from .pof_autogen import * 
from .intfactor import * 
from morebs2.numerical_generator import prg_seqsort_ties,prg_partition_for_sz

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

    def add_to_acc_queue(self,q): 
        self.acc_queue.extend(q)  
    
    def update_acc_queue(self,s): 
        self.acc_queue = list(set(self.acc_queue) - set(s))

    #----------------------- children-access functions 
    """
    def add_children(self,idn_seq,entryf_seq): 
        assert len(idn_seq) == len(entryf_seq)
        for (i,x) in enumerate(idn_seq): 
            q = IDecNode(x,entryf_seq[i],None,self.rd+1)
            self.children.append(q)
    """

    def add_children_nodes(self,tnx): 
        for tn in tnx: 
            assert type(tn) == IDecNode
            self.children.append(tn)

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
        return self.travf.apply(v) 

"""
represents a factor-based or polynomial-based boolean classifier 
"""
class IDecTrFunc: 

    def __init__(self,conditional_list,conditional_type,switched_indices=None):  
        assert conditional_type in {"poly","factor"}
        if conditional_type == "poly": 
            for x in conditional_list: 
                assert type(x) in {PolyOutputFitterVar1,PolyOutputFitterVar2}
        else: 
            for x in conditional_list: 
                assert type(x) in {int,np.int32,np.int64} 


        if type(switched_indices) == type(None): 
            switched_indices = [0] * len(conditional_list) 
        assert set(switched_indices).issubset({0,1})
        assert len(switched_indices) == len(conditional_list)

        self.cl = conditional_list
        self.ct = conditional_type
        self.si = switched_indices
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
        stat = x / self.cl[i] == x // self.cl[i]
        if self.si[i]: 
            return not stat 
        return stat 

    def poly_classify(self,x,i): 
        assert self.ct == "poly" 
        stat = False 
        if type(self.cl[i]) == PolyOutputFitterVar1:
            try: 
                x1,x2 = self.cl[i].apply(x), self.cl[i].c  
                stat = x1 == x2 
            except: 
                stat = False 
        else: 
            try: 
                x1,x2 = self.cl[i].apply(x), self.cl[i].apply(self.cl[i].x1)  
                stat = x1 == x2 
            except:
                stat = False 

        if self.si[i]: 
            return not stat 
        return stat 

    def switch_conditional(self,ci): 
        self.si[ci] = (self.si[ci] + 1) % 2 

    def one_switch(self,x): 
        assert x in {0,1} 
        self.si = [x] * len(self.si)

"""
every <IDecNode>'s `travf` variable is an instance of this, used 
for directing input values towards next nodes. A multi-classifier 
built from an ordered sequence of <IDecTrFunc> instances. 
"""
class IDecNodeTravFunc:

    def __init__(self):
        self.bclassif = [] 
        self.bclassif_nextnode = []
        self.cat_samples = []  
        self.default_class = None
        self.default_class_samples = None  
        return

    def __str__(self): 
        s = "pos-labels: " + str(len(self.bclassif)) + "\n"
        for (i,cs) in enumerate(self.cat_samples): 
            if type(cs) == type(None): 
                l = 0 
            else: 
                l = len(cs) 
            s += "label: {} \tsize: {}".format(self.bclassif_nextnode[i],l) + "\n"
        return s 

    def set_default_class(self,dc,samples):
        assert type(dc) in {int,np.int32,np.int64}
        self.default_class = dc 
        self.default_class_samples = samples 

    def add_bclassif_nextnode_pair(self,bc,nn,cat_sample=None): 
        assert type(bc) == IDecTrFunc
        assert nn not in self.bclassif_nextnode
        assert type(cat_sample) in {list,set,type(None)}
        self.bclassif.append(bc)
        self.bclassif_nextnode.append(nn)
        self.cat_samples.append(cat_sample)

    def apply(self,x): 
        for (i,c) in enumerate(self.bclassif): 
            if c.bclassify(x): 
                return self.bclassif_nextnode[i] 
        return self.default_class

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
        self.node_ctr = 0 
        self.node_cache = [] 
        self.root = None 

    def factor_preproc(self):
        self.isfso = ISFactorSetOps(deepcopy(self.intseq.l),\
            int_limit=DEFAULT_INT_MAX_THRESHOLD)
        self.isfso.factor_count_()
        return 

    #----------------------- root initialization to satisfy depth or
    #----------------------- leaf requirement 

    def init_root(self): 
        tn = IDecNode(self.node_ctr,None,None,0) 
        tn.add_to_acc_queue(deepcopy(self.intseq.l))
        self.root = tn 
        self.node_cache.append(tn) 
        self.node_ctr += 1 
        return

    def init_lf(self): 
        return -1 

    def init_df(self):
        return -1 

    def split_node(self,node): 
        return -1

    #---------------------- splitter for depth 

    """
    satisfies the depth requirement `d`, if not None, by 
    declaring nodes that form a path of length `d` starting 
    from `tn` (typically the root node). The path constitutes 
    a 1-to-1 mapping of the `d` elements to the `d` nodes 
    (every element satisfies exactly one node). 
    """
    def split__depthreq(self,tn,split_type): 
        assert split_type in {"poly","factor"}

        S = deepcopy(tn.acc_queue) 

        # declare the first split
        if split_type == "poly": 
            classif,siblings = self.poly_subset_classifier(S,self.d)
            idntf = IDecNodeTravFunc()
            idntf.add_bclassif_nextnode_pair(classif,self.node_ctr,siblings)
            self.node_ctr += 1 
            tn.set_travf(idntf) 
            tn.update_acc_queue(siblings)
        else: 
            prt0 = [self.d,len(S) - self.d] 
            travf = self.factor_split_travf(S,prt0,last_subset_isneg=True)
            tn.set_travf(travf)
            tn.update_acc_queue(travf.cat_samples[0])

        cx = IDecNode(tn.travf.bclassif_nextnode[0],None,None,tn.rd+1) 
        cx.add_to_acc_queue(tn.travf.cat_samples[0])
        tn.add_children_nodes([cx]) 

        S = cx.acc_queue
        stat = len(cx.acc_queue) > 1
        while stat: 
            S = cx.acc_queue
            if split_type == "poly":
                travf = self.poly_one_classify(S) 
                aqueue = travf.cat_samples[0]
            else: 
                prt0 = [1,len(S) - 1]  
                travf = self.factor_split_travf(S,prt0,last_subset_isneg=True)

                # all elements except for the fitted one can pass 
                travf.bclassif[0].switch_conditional(0) 
                aqueue = list(set(S) - travf.cat_samples[0])

            cx2 = IDecNode(travf.bclassif_nextnode[0],None,None,cx.rd+1)
            cx2.add_to_acc_queue(aqueue)
            cx.add_children_nodes([cx2]) 
            cx.set_travf(travf) 

            cx = cx2
            stat = len(cx2.acc_queue) >= 2 

        if len(tn.acc_queue) > 0: 
            self.node_cache.append(tn) 

    #---------------------- factor-splitter functions 
    
    def factor_split_travf(self,S,partition,last_subset_isneg:bool=False,\
        set_default_class:bool=False): 
        fspart = self.factor_split__partitioned(S,partition,\
            last_subset_isneg)
        idntf = IDecNodeTravFunc()

        S_ = set(S) 
        for (k,fx) in fspart.items(): 
            classif = self.partition_subset_to_factor_classifier(fx)
            samples = set() 
            for fx_ in fx: samples |= fx_[1] 
            next_node = self.node_ctr 
            self.node_ctr += 1
            idntf.add_bclassif_nextnode_pair(classif,next_node,samples)
            S_ -= samples 

        if set_default_class: 
            idntf.set_default_class(self.node_ctr,list(S_))
            self.node_ctr += 1 

        return idntf 

    def partition_subset_to_factor_classifier(self,ps): 
        conditional_list = [ps_[0] for ps_ in ps] 
        return IDecTrFunc(conditional_list,"factor")

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
            fx,S = self.factors_for_subset(S,p) 
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
            S = S - keys 
            self.update_sorted_factors(keys)
            #self.isfso.remove_seq_elements(keys) 
            f.append((factor[0],keys))  
            sz -= len(keys) 
        return f,S

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

    #---------------------- poly-splitter functions 
    
    # TODO: test 
    def poly_subset_classifier(self,S,class_size:int):
        assert len(S) >= class_size
        assert class_size >= 2
        assert len(np.unique(S)) >= 2 

        # choose two elements from S to pair up
        j = self.prg() % len(S) 
        x1 = S.pop(j)
        x2 = None 
        while True: 
            j = self.prg() % len(S) 
            if S[j] != x1: 
                x2 = S.pop(j)
                break 

        S2 = set([x1,x2]) 

        # solve the polynomial 
        pofgen = POFV2ConditionAutoGen(self.prg) 

        pofv2 = pofgen.integerpair_op(x1,x2,\
            sibling_range=[1,2],coeff_range=DEFAULT_COEFF_RANGE,\
            power_range = DEFAULT_POWER_RANGE,deepcopy_prng=False)
        pofv2 = next(pofv2) 

        while len(S2) < class_size:
            j = self.prg() % len(S) 
            sx = S.pop(j) 
            S2 |= {sx}

        S2 -= {x1,x2}
        pofv1_siblings = pofgen.POFV2_to_POFV1_siblings(pofv2,S2) 

        conditional_list = [pofv2] + pofv1_siblings
        return IDecTrFunc(conditional_list,"poly"), S2 | {x1,x2} 

    def poly_one_classify(self,S):

        # choose a value from S 
        j = self.prg() % len(S)
        x = S.pop(j)

        # choose a value in the range of [1,10**6]
        q = modulo_in_range(self.prg(),DEFAULT_COEFF_RANGE) 
        q = 1 if q == 0 else q 

        # choose an exponent
        n = modulo_in_range(self.prg(),DEFAULT_POWER_RANGE)

        pofv1 = PolyOutputFitterVar1(n,x,q,self.prg,default_sizemod=False)
        pofv1.solve()

        conditional_list = [pofv1] 
        idtf = IDecTrFunc(conditional_list,"poly")
        idtf.one_switch(1)

        idntf = IDecNodeTravFunc()
        idntf.add_bclassif_nextnode_pair(idtf,self.node_ctr,set(S)) 
        return idntf 