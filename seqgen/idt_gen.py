from intigers.mod_prng import *
from intigers.idt_proc import * 
from intigers.extraneous import * 
from morebs2.numerical_generator import prg_seqsort,prg_choose_n

DEFAULT_FOREST_NEWSEQ_NUMCENTERS = [2,15]
DEFAULT_FOREST_NEWSEQ_MULTRANGE = [-20,20]

class IDecForest: 

    def __init__(self,s,prng_outputter,cache_size,reprod_rate_range,\
        tree_cap:int,prg,prg2=None,verbose:bool = True):
        assert type(s) == IntSeq
        assert len(s) >= 5 
        assert type(prng_outputter) == ModPRNGOutputter
        assert reprod_rate_range[0] < reprod_rate_range[1] 
        assert reprod_rate_range[0] > 0 
        assert type(reprod_rate_range[0]) == type(reprod_rate_range[1])
        assert type(reprod_rate_range[0]) == int 
        assert tree_cap > 0 and type(tree_cap) == int 

        self.ST = []
        self.S = s
        
        self.l = len(self.S)
        self.default_length_range = \
            IDecForest.default_IDecForest_sequence_range(self.l)
        self.default_integer_range = \
            IDecForest.default_IDecForest_integer_range(min(self.S),max(self.S))
        self.T = None 
        self.mpo = prng_outputter
        # ? 
        self.queue = []
        self.cache_size = cache_size
        self.reprod_rate_range = reprod_rate_range
        self.tree_cap = tree_cap 
        self.prg = prg 
        self.prg2 = prg2
        self.verbose = verbose 
        self.next_reprod = None 
        self.reprod_ctr = None
        self.reset_reprod_ctr()

    """
    main method: test 
    """
    def __next__(self): 

        if len(self.queue) == 0: 
            # case: initialize
            if len(self.ST) == 0: 
                self.one_tree() 

            # choose a tree
            q = self.prg() % len(self.ST)
            T = self.ST[q][1] 
            if self.verbose: print("choose tree")

            # choose a sequence to process
            seqtype = self.prg() % 3
            S = self.one_new_IntSeq(seqtype)
            
            if self.verbose: print("new seq: ",S) 

            # choose a process type
            proctype = self.prg() % 4

            # add output sequence to queue
            q = self.process_seq_at_tree(T,S,proctype)
            self.queue.extend(q) 

        x = self.queue.pop(0)
        self.reprod_ctr += 1

        # case: reproduce new tree
        if self.reprod_ctr >= self.next_reprod:
            seqtype = self.prg() % 3
            self.S = self.one_new_IntSeq(seqtype)
            if self.verbose: 
                print("reproducing @: ",self.reprod_ctr)
                print("seq: ",self.S)
            self.one_tree()
            self.reset_reprod_ctr()

        return x

    def reset_reprod_ctr(self): 
        self.reprod_ctr = 0
        # reset `next_reprod`
        if len(self.ST) >= self.tree_cap:
            self.next_reprod = float('inf') 
        else: 
            self.next_reprod = modulo_in_range(\
                self.prg(),self.reprod_rate_range)

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
        
        cnvrt = IntSeq2Tree(I,l,d,prgx,verbose=False)
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
            q = intlist_no_dups_no_zero_abs(q,px)
            q = IntSeq(q) 
        else: 
            # choose a sequence 
            assert len(self.ST) > 0
            i = self.prg() % len(self.ST)
            S = self.ST[i][0] 
            q = IDecForest.derive_new_IntSeq__type_fm(S,px,self.default_length_range)
            q = IntSeq(intlist_no_dups_no_zero_abs(q.l,px))
            
            # case: add noise
            if gentype == 2:
                l = len(q)
                v = [modulo_in_range(px(),self.default_integer_range) \
                    for _ in range(len(q))]
                v = np.array(v,dtype=np.int32)
                q = q + v
                v2 = [modulo_in_range(v_,self.default_integer_range) \
                    for v_ in q]
                v2 = intlist_no_dups_no_zero_abs(v2,px)
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

        iseq = intlist_no_dups_no_zero_abs(iseq,prg)
        return IntSeq(iseq)

    #------------------------ output generation

    def process_seq_at_tree(self,T,S,proc_type):
        assert proc_type in {0,1,2,3}

        if self.verbose: print("proc type: ",proc_type)
        if proc_type == 1:
            # make new sequence 
            seqtype = self.prg() % 4
            intseq = self.one_new_IntSeq(seqtype)

            if len(intseq) > len(S):
                S1,S2 = S,intseq
            else: 
                S1,S2 = intseq,S 
            return self.process_seq_at_tree__iso_sequential(T,S1,S2)

        elif proc_type == 0: 
            return self.process_seq_at_tree__sequential(T,S)
        elif proc_type == 2: 
            return self.process_seq_at_tree__inflow(T,S)
        else:
            return self.process_seq_at_tree__splat(T,S)

    # choose seq for tree 
    def process_seq_at_tree__sequential(self,T,S):
        assert type(S) == IntSeq
        itp = IDTProc(T)

        qs = []
        for s in S: 
            _,vs = itp.process_value(s,True)
            qs.extend(vs)
        return qs 

    def process_seq_at_tree__iso_sequential(self,T,S,S2):

        assert type(S) == IntSeq
        itp = IDTProc(T)

        qs = []
        for (i,s2) in enumerate(S2): 
            s1 = S[i % len(S)]
            qs_ = itp.iso_output(s1,s2)
            qs.extend(qs_)
        return qs 

    # TODO: test 
    def process_seq_at_tree__inflow(self,T,S):

        # form the partition
        prt = IDecForest.default_IDecForest_partitioning(len(S),\
            self.prg,self.prg2)

        # shuffle the sequence
        i2 = list(S.l)
        i2 = prg_seqsort(i2,self.prg)  

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
            qr,stat = self.process_one_partition_set__inflow(itp,\
                set(st),num_rounds=dt)
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
            qrx = next(itp)
            if len(qrx) == 0: 
                stat = True
                break
            qr_ = IDecForest.default_IDecForest_shuffle_dictoutput(qrx,self.prg)
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
        SD = IDecForest.default_IDecForest_splatdict(nodelist,intlist,self.prg)

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