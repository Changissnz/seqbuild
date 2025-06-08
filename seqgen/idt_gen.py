from intigers.mod_prng import *
from intigers.idt_proc import * 
from intigers.extraneous import * 

DEFAULT_FOREST_NEWSEQ_NUMCENTERS = [2,15]
DEFAULT_FOREST_NEWSEQ_MULTRANGE = [-20,20]

class IDecForest: 

    def __init__(self,s,prng_outputter,cache_size,prg,prg2=None):
        assert type(s) == IntSeq
        assert len(s) >= 5 
        assert type(prng_outputter) == ModPRNGOutputter
        self.ST = []
        self.S = s
        
        self.l = len(self.S)
        self.default_length_range = \
            IDecForest.default_IDecForest_sequence_range(self.l)
        self.default_integer_range = \
            IDecForest.default_IDecForest_integer_range(min(self.S),max(self.S))
        self.T = None 
        self.mpo = prng_outputter
        self.cache_size = cache_size
        self.prg = prg 
        self.prg2 = prg2

    """
    outputs n values
    """
    def output_n(self,n):  
        return -1

    #---------------------- tree and integer sequence creation

    def one_tree(self): 
        self.one_tree_()
        self.ST.append((self.S,self.T))

    def one_tree_(self): 
        # fetch arguments for <IntSeq2Tree> conversion 
        I = self.S
        
        l,d = None,None 
        rnge = [2,ceil(len(I) / 2)] 
        if rnge[0] == rnge[1]: rnge[1] += 1 

        q = modulo_in_range(self.prg(),rnge)
        if self.prg() % 2: 
            l = q
        else: 
            d = q     

        prgx = self.prg2 if type(self.prg2) != type(None) else self.prg 
        cnvrt = IntSeq2Tree(I,l,q,prgx)
        cnvrt.convert()

        # set the `entryf` functions for the new tree 
        self.T = cnvrt.root 
        TNode.dfs(self.T,False,False,True,set_attr=('entryf',self.mpo.__next__))
        return self.T

    """
    outputs one new <IntSeq> by one of three generative types: 
    - 0 := values from `prg2` OR `prg`
    - 1 := factor-multiple derivation of existing sequence
    - 2 := factor-multiple derivation of existing sequence with noise from 
           `prg2` or `prg`. 
    """
    def one_new_IntSeq(self,gentype:int):
        assert gentype in {0,1,2} 

        # from prg2 OR prg 
        px = self.prg2 if type(self.prg2) != type(None) \
            else self.prg

        q = None 
        if gentype == 0: 
            l2 = modulo_in_range(px(),self.default_length_range)
            q = []
            for _ in range(l2):
                n2 = modulo_in_range(px(),self.default_integer_range)
                q.append(n2) 
            q = intlist_no_dups_no_zero(q)
            q = IntSeq(q) 
        else: 
            # choose a sequence 
            assert len(self.ST) > 0
            i = self.prg() % len(self.ST)
            S = self.ST[i][0] 
            q = IDecForest.derive_new_IntSeq__type_fm(S,px,self.default_length_range)

            # case: add noise
            if gentype == 2:
                l = len(q)
                v = [modulo_in_range(px(),self.default_integer_range) \
                    for _ in range(len(q))]
                q = q + v
                v2 = [modulo_in_range(v_,self.default_integer_range) \
                    for _ in range(len(q))]
                q = IntSeq(v2) 
        return q 

    """
    """
    @staticmethod
    def derive_new_IntSeq__type_fm(S,prg,length_range):
        assert type(S) == IntSeq

        # set number of centers
        ncr = (DEFAULT_FOREST_NEWSEQ_NUMCENTERS[0],\
            min([len(S),DEFAULT_FOREST_NEWSEQ_NUMCENTERS[1]]))
        c_ = modulo_in_range(prg(),ncr)

        # fetch the centers, each a factor
        isfso = ISFactorSetOps(S.l,int_limit=NPINT32_MAX)
        isfso.factor_count_()
        q = isfso.dsort(pkeys=None)

        cx = []
        while c_ > 0:
            i = prg() % len(q)
            cx.append(q.pop(i)[0])
            c_ -= 1 

        # generate the sequence 
        n = modulo_in_range(prg(),length_range)

        rx = modulo_in_range(prg(),[0,1001])
        rx = rx / 1000.0 

        iseq,_ = prg__integer_seq__mult(n,cx,rx,None,\
            DEFAULT_FOREST_NEWSEQ_MULTRANGE,prg,\
            num_attempts_per_nc=150)

        iseq = intlist_no_dups_no_zero(iseq)
        return IntSeq(iseq)

    #------------------------ output generation

    # choose seq for tree 
    def process_seq_at_tree__sequential(self,T,S):
        assert type(S) == IntSeq
        itp = IDTProc(T)

        qs = []
        for s in S: 
            _,vs = itp.process_value(s)
            qs.extend(vs)
        return qs 

    def process_seq_at_tree__iso_sequential(self,T,S,S2):

        assert type(S) == IntSeq
        itp = IDTProc(T)

        qs = []
        for (i,s2) in enumerate(S2): 
            s1 = S[i % len(S)]
            qs_ = itp.iso_output(v1,s2)
            qs.extend(qs_)
        return qs 

    # TODO: test 
    def process_seq_at_tree__inflow(self,T,S):

        # form the partition
        prt = IDecForest.default_IDecForest_partitioning(len(S),\
            self.prg,self.prg2)

        # shuffle the sequence
        i2 = list(S.l)
        i2 = prg_seqsort(i2,px) 

        # get the max frequency for calling `inflow_set` 
        _,depth,_ = TNode.dfs(T,display=False,\
            collect=True,reset_index=True,set_attr=None,\
            fetch_set=set())
        depth_range = [1,depth+1]
        if depth_range[0] == depth_range[1]: 
            depth_range[1] += 1 

        # iterate through the partitions
        itp = IDTProc(T)

        s = 0 
        qrs = [] 
        for p in prt:
            st = i2[s:s+p]
            dt = modulo_in_range(self.prg(),depth_range)
            qr,stat = self.process_one_partition_set__inflow(itp,st,num_rounds=dt)
            qrs.extend(qr) 
            s = s + p

        while len(itp.flow_queue) > 0: 
            qr = next(itp)
            if len(qr) == 0: break
            qr_ = IDecForest.default_IDecForest_shuffle_dictoutput(qr,self.prg)
            qrs.extend(qr_)
        return qrs 

    def process_one_partition_set__inflow(self,itp,s,num_rounds=float('inf')): 
        assert num_rounds > 0 

        c = 0
        itp.inflow_set(s) 
        qr = [] 
        stat = False
        while c < num_rounds:
            qr = next(itp)
            if len(qr) == 0: 
                stat = True
                break
            qr_ = IDecForest.default_IDecForest_shuffle_dictoutput(qr,self.prg)
            qr.extend(qr_)
            c += 1 

        return qr,stat 

    # TODO: test 
    def process_seq_at_tree__splat(self,T,S):
        D,_,_ = TNode.dfs(T,False,True,True,set_attr=None)

        # choose between 25 and 75 percent of the samples from S
        q = modulo_in_range(self.prg(),[0,1001]) / 1000.0 
        n = int(ceil(len(S) / 2.0))  

        nodelist = sorted(list(D.keys()))
        intlist = prg_choose_n(list(S.l),n,self.prg,is_unique_picker=True)
        SD = default_IDecForest_splatdict(nodelist,intlist,self.prg)

        # splat now 
        itp = IDTProc(T) 
        return itp.splat_process(SD)

    #------------------------- default methods for specifying ranges 
    #------------------------- of sequence lengths and integer values, 
    #------------------------- amongst other tasks. 

    """
    Default range for an integer sequence generated by a PRNG. 
    This default range is used in the case of generative type #0 
    (see method `one_new_IntSeq`). 
    """
    @staticmethod 
    def default_IDecForest_sequence_range(start_sequence_length):
        minnie,maxie = int(ceil(start_sequence_length / 2)),\
            start_sequence_length * 2
        return (minnie,maxie) 

    @staticmethod
    def default_IDecForest_integer_range(smin,smax):
        minnie,maxie = int(ceil(smin / 2)),\
            smax * 2
        return (minnie,maxie)  

    @staticmethod
    def default_IDecForest_partitioning(l,prg,prg2):

        mx = int(ceil(l / 2))
        assert mx >= 2
        if mx == 2: mx += 1
        num_sets = modulo_in_range(prg(),[2,mx])

        px = prg2 if type(prg2) != type(None) \
            else prg 
        var = modulo_in_range(prg(),[0,1001])
        var = var / 1000.0

        prt = prg_partition_for_sz(l,num_sets,px,var)
        return prt 

    """
    """
    @staticmethod
    def default_IDecForest_shuffle_dictoutput(D,prg):
        vs = []
        ks = sorted(D.keys())

        for k in ks:
            vs.extend(D[k])
        return prg_seqsort(vs,prg)

    @staticmethod
    def default_IDecForest_splatdict(nodelist,intlist,prg):
        D = defaultdict(list)

        for i in intlist: 
            q = prg() % len(nodelist) 
            n = nodelist[q] 
            D[n].append(i)
        return D 