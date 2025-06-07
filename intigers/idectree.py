from .pof_autogen import * 
from .intfactor import * 
from morebs2.numerical_generator import prg_seqsort_ties,prg_partition_for_sz
from morebs2.g2tdecomp import TNode 

"""
`p` is a partition that does not contain 0. Fixes any 
element of `p` that is less than 2 by transferring 1 
over from an element of value greater than 2. 

:return: p' if p has any elements below 2 o.w. p, success status 
        of fix. 
:rtype: list|np.ndarray,bool. 
"""
# NOTE: method used for partitions valid to the requirements 
#       of polynomial-splitting. 
def partition_fix__subset_is_minsize_2(p,prg): 
    p_ = np.array(p) 
    indices = np.where(p_ < 2)[0]
    if len(indices) == 0: return p_,True 

    other_indices = np.where(p_ > 2)[0]
    stat = True 
    for i in indices: 
        if len(other_indices) == 0: 
            stat = False 
            break 

        diff = 2 - p_[i]
        j2 = prg() % len(other_indices) 
        j = other_indices[j2] 
        p_[j] -= 1 
        p_[i] += 1 

        if p_[j] == 2:
            other_indices = np.delete(other_indices,j2)
    return p_,stat 

#------------------------------------------------------------------
"""
Node structure used as unit for a decision tree constructed 
from an <IntSeq>. Structure is associated with functions, 
an input (operation) function `entryf` for an incoming value 
transmitted by a parent and an output (classifier) function `travf` 
to map the incoming value to the next node. 
"""
class IDecNode(TNode): 

    def __init__(self,idn:int,entryf,travf,root_distance:int): 
        assert type(idn) == int  
        is_root = True if root_distance == 0 else False
        super().__init__(idn,False,is_root,root_distance)

        self.entryf = entryf
        self.travf = travf 
        self.children = [] 
        self.reachability_key = None 
        self.con_key = set() 
        self.acc_queue = []
        return

    #----------------------- var-setter functions 

    def full_clear(self):
        qx = [self]

        while len(qx) > 0:
            qx2 = []
            for x in qx:
                x.clear()
                qx2.extend(x.children)
            qx = qx2

    def clear(self): 
        self.acc_queue.clear() 
        if type(self.travf) != type(None): 
            self.travf.clear() 

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

    def add_children(self,tnx): 
        for tn in tnx: 
            assert type(tn) == IDecNode
            self.children.append(tn)

    def fetch_conn(self,idn): 
        i = self.index_of_child(idn) 
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

    def __add__(self,x):
        assert type(x) == IDecNodeTravFunc

        assert set(self.bclassif_nextnode).intersection(set(\
            x.bclassif_nextnode)) == set() 
        q = deepcopy(self)
        q.bclassif.extend(x.bclassif)
        q.bclassif_nextnode.extend(x.bclassif_nextnode)
        q.cat_samples.extend(x.cat_samples)
        if type(x.default_class) != type(None): 
            q.default_class = x.default_class
            q.default_class_samples = x.default_class_samples
        return q  


    def __str__(self): 
        s = "pos-labels: " + str(len(self.bclassif)) + "\n"
        for (i,cs) in enumerate(self.cat_samples): 
            if type(cs) == type(None): 
                l = 0 
            else: 
                l = len(cs) 
            s += "stype: {} \tlabel: {} \tsize: {}".format(\
                self.bclassif[i].ct,self.bclassif_nextnode[i],l) + "\n"
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

    def clear(self):
        self.cat_samples.clear()  


# NOTE: some code in this class may need to be refactored. 
"""
converts an <IntSeq> instance to a tree (directed graph) T. 
Specifically, the tree is single root with single-parented 
children.  

T satisfies a leaf requirement `l` XOR  a depth requirement `d`. 

T starts as a single node that takes as input the the integer 
sequence `intseq`. A process calculates a classifier to split 
the sequence according to an arbitrary value valid in accordance 
with the size of the integer sequence. The splitting of a sequence 
at a node results in disjoint subsets of the sequence getting routed 
to different destinations, namely at the node (value is no longer 
considered in splitting process) or to a child node of the node. 
The classifier is an instance of <IDecNodeTravFunc>, which is comprised 
of an ordered sequence of boolean classifier instances <IDecTrFunc>. 
Each of these boolean classifiers C checks an input value v to see if 
it satisfies any of its conditionals c_i: 
        C(v) = 1 if there exists a c_i in C such that c_i(v).  
For a boolean classifier C, all conditionals are of one type `factor` 
or `poly`. The type `factor` splits sequences according to the sequence 
elements that are multiples of the corresponding factor. For the `poly` 
type, at least two integers are required for the initial fitting via 
a <PolyOutputFitterVar2> instance. Once these two integers are fitted 
for the first conditional, any remaining number of integers can be added 
to the <IDecTrFunc> by fitting them with a <PolyOutputFitterVar1> to the 
output value from the <PolyOutputFitterVar2>. To split, the `factor` type 
requires at least two values, and the `poly` type requires at least four 
values. A boolean classifier corresponds to a node identifier, one child 
of the node the classifier belongs to. 

An <IDecNodeTravFunc> iterates through its sequence of boolean classifiers 
in order to classify an input value v. If a boolean classifier C outputs 1, 
then v is of the category C corresponds to. If all boolean classifiers output 
0, exactly one of two things happens:
- the value v ceases travel and its final destination is the node it is at, 
- the value v travels to the node of identifier `default_class`, the class 
  variable of <IDecNodeTravFunc>, if `default_class` is not None. 

The conversion process starts with satisfying either the leaf or depth 
requirement. A pseudo-random number (integer) generator `prg` is used for making 
certain decisions in this first part. For the remaining elements of the sequence, 
the `prg` outputs values that decide attributes on splitting, such as: 
- partition for the elements
- type of boolean classifier

For all nodes except for the root, their <IDecNodeTravFunc> function will 
always have instances of <IDecTrFunc> that are all of exactly `factor` or 
`poly`, due to the programming. The root `travf` may have <IDecTrunc> functions 
of differing types. 

NOTE: some deficiencies of this splitting algorithm, due to the use of factors 
      and polynomial-based equations, are 
- duplicate integers can never be split against each other. 
- the value 0 cannot be used for splitting. 
- (infrastructure-dependent) numbers larger than ?3200? will break the program. 
"""
class IntSeq2Tree: 

    def __init__(self,intseq,l:int,d:int,prg,verbose:bool=False): 
        assert type(intseq) == IntSeq
        assert type(l) == type(None) or type(d) == type(None)
        self.intseq = intseq 
        self.l = l 
        self.d = d
        self.verbose=verbose
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
        if self.verbose:
            print("finished factor count.")
        return 

    """
    main method 
    """
    def convert(self): 
        self.init_root()
        if self.leaf_first: 
            self.init_lf()
        else: 
            self.init_df()
        if self.verbose: 
            print("splitting the rest.")
            print("---/---/---/---/")

        while len(self.node_cache) > 0: 
            tn = self.node_cache.pop(0)
            self.split_node(tn,partition=None)
            if self.verbose: 
                print('---/---/---/---/')

    #----------------------- root initialization to satisfy depth or
    #----------------------- leaf requirement 

    def init_root(self): 
        tn = IDecNode(self.node_ctr,None,None,0) 
        tn.add_to_acc_queue(deepcopy(self.intseq.l))
        self.root = tn 
        self.node_cache.append(tn) 
        self.node_ctr += 1 
        if self.verbose:
            print("- init root.")
        return

    def init_lf(self): 
        tn = self.node_cache.pop(0) 

        num_sets = self.l 
        partition = self.partition_for_node(tn,num_sets=num_sets)
        self.split_node(tn,partition=partition) 
        return

    def init_df(self):
        tn = self.node_cache.pop(0) 
        f = "factor" if self.prg() % 2 else \
            "poly"

        self.split__depthreq(tn,f)
        return

    #---------------------- main splitting function

    def split_node(self,node,partition=None):
        if len(node.acc_queue) < 2:
            return 

        # get the split type 
        is_factor = False 
        if self.prg() % 2: 
            is_factor = True 

        if len(node.acc_queue) < 4:
            is_factor = True 
        
        # get the partition 
        is_min2 = True if not is_factor else False 
        if type(partition) == type(None): 
            partition = self.partition_for_node(node,is_min2=is_min2) 
        else: 
            assert sum(partition) == len(node.acc_queue)
            if is_min2: 
                partition,stat = partition_fix__subset_is_minsize_2(partition,self.prg)
                assert stat 

        if self.verbose: 
            stype = "factor" if is_factor else "poly"
            print("\t- split type {}, partition: {}".format(\
                stype,partition))

        if is_factor: 
            travf = self.factor_split_travf(deepcopy(node.acc_queue),\
                partition,last_subset_isneg=False,set_default_class=False)
        else: 
            travf = self.poly_split__travf(deepcopy(node.acc_queue),partition,\
                last_subset_isneg=False,set_default_class=False)

        self.set_travf_for_node(node,travf) 

    """
    sets the output function `travf` and the corresponding children 
    nodes for `node`. 
    """
    def set_travf_for_node(self,node,travf): 

        # declare the children 
        tns = [] 
        for cl in travf.bclassif_nextnode: 
            tn = IDecNode(cl,None,None,node.rdistance+1) 
            tns.append(tn)
        
        # check for default class
        if type(travf.default_class) != type(None): 
            tn = IDecNode(travf.default_class,None,None,node.rdistance+1) 
            tns.append(tn) 
        node.add_children(tns) 

        # classify all samples
        D = defaultdict(list)
        while len(node.acc_queue) > 0: 
            x = node.acc_queue.pop(0) 
            ci = travf.apply(x)
            if type(ci) == type(None): 
                continue 
            D[ci].append(x) 

        # add samples to children pertaining to 
        # their category 
        for k,v in D.items(): 
            c = node.fetch_conn(k) 
            c.add_to_acc_queue(v) 
            if len(c.acc_queue) > 1: 
                self.node_cache.append(c)

        if type(node.travf) != type(None): 
            travf = node.travf + travf 

        node.set_travf(travf) 
        node.clear() 

    def partition_for_node(self,node,num_sets=None,is_min2=False): 
        # choose the number of sets + variance for a partition
        if self.verbose: 
            print("partitioning {} elements".format(\
                len(node.acc_queue)))
            i = 0
            while True:
                aq0 = node.acc_queue[i * 5: (i+1) * 5]
                print("\t {}".format(aq0))
                i += 1 
                if i * 5 > len(node.acc_queue): break 
            
        assert len(node.acc_queue) > 1 

        if type(num_sets) != type(None): 
            assert num_sets <= len(node.acc_queue)
        else: 
            max_sets = ceil(len(node.acc_queue) / 2)# + 1 
            max_sets = max([max_sets,3])
            num_sets = modulo_in_range(self.prg(),[2,max_sets])

        variance = modulo_in_range(self.prg(),[0,1001]) / 1000.0 
        partition = prg_partition_for_sz(len(node.acc_queue),num_sets,\
            self.prg,variance)

        # NOTE: caution. 
        if is_min2: 
            partition,stat = partition_fix__subset_is_minsize_2(partition,self.prg) 
            if not stat: raise ValueError("something wrong.")

        return partition

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

        if self.verbose:
            print("start split for depth-req")

        # declare the first split
        if split_type == "poly": 
            classif,siblings = self.poly_subset_bclassifier(S,self.d)
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

        cx = IDecNode(tn.travf.bclassif_nextnode[0],None,None,tn.rdistance+1) 
        cx.add_to_acc_queue(tn.travf.cat_samples[0])
        tn.add_children([cx]) 

        if self.verbose: print("\t- first split of ({}) done.".format(split_type))
    
        S = cx.acc_queue
        stat = len(cx.acc_queue) > 1
        while stat:             
            S = cx.acc_queue
            if self.verbose: print("\t* remaining: ",len(S))
            if split_type == "poly":
                travf = self.poly_one_classify(S) 
                aqueue = travf.cat_samples[0]
            else: 
                prt0 = [1,len(S) - 1]  
                travf = self.factor_split_travf(S,prt0,last_subset_isneg=True)

                # all elements except for the fitted one can pass 
                travf.bclassif[0].switch_conditional(0) 
                aqueue = list(set(S) - travf.cat_samples[0])

            cx2 = IDecNode(travf.bclassif_nextnode[0],None,None,cx.rdistance+1)
            cx2.add_to_acc_queue(aqueue)
            cx.add_children([cx2]) 
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

    def poly_split__travf(self,S,partition,last_subset_isneg:bool=False,\
        set_default_class:bool=False):

        l = len(partition) if not last_subset_isneg \
            else len(partition) - 1

        travf = IDecNodeTravFunc() 
        for i in range(l): 
            p = partition[i] 
            bfunc,q = self.poly_subset_bclassifier(S,p)
            S = list(set(S) - q)
            travf.add_bclassif_nextnode_pair(bfunc,self.node_ctr)
            self.node_ctr += 1

        if set_default_class: 
            travf.set_default_class(self.node_ctr,None) 
            self.node_ctr += 1 

        return travf 
    
    def poly_subset_bclassifier(self,S,class_size:int):
        assert len(S) >= class_size
        assert class_size >= 2
        assert len(np.unique(S)) >= 2 

        if self.verbose: 
            print("\t\tpoly classifiers for set, size {},subset {}".\
                format(len(S),class_size))

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

        # make siblings for the polynomial 
        pofv1_siblings,solvestat = pofgen.POFV2_to_POFV1_siblings(pofv2,S2) 
        if self.verbose: print("poly-split number of siblings: ",len(S2) + 2)

        # ensure all siblings solved constant c 
        assert set(solvestat) == {True} or len(S2) == 0 

        # make the `travf` function 
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
        self.node_ctr += 1 
        return idntf 