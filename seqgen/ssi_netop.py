from .ssi_load import * 
from morebs2.numerical_generator import prg_partition_for_sz__n_rounds
from morebs2.g2tdecomp import * 

"""
(S)earch (S)pace (I)terator Net Operator. 

A pseudo-random number generator that relies in part on structure provided by 
a search space iterator. 

        *Instantiation Variables* 
Set `io_prg` to a PRNG to add noise to pre-output values that are responsible 
for the output values from <SSINetOp>. 

Set `shuffle_dist` to 1 (active) for the <SSINetOp> to shuffle the numbers in 
its reference queue Q, every time before using the shuffled numbers to update 
the pertinent <SSINetNode__TypeLCGNet>s. 

The reference queue Q is `tmp_queue` (contains numbers from all structure types 
`lcg`,`mdr*`,`optri`) if `lcg_input_only` is set to 0. Otherwise, Q is `lcg_queue`. 

        *Qualities* 

In a typical instantiation via method<SSINetOp.one_instance>, <SSIBatchLoader__TypeLCGNet> 
is used to generate at least five <SSINetNode__TypeLCGNet> instances. These instances contain 
at least one of each of the following structure types: "lcg","mdr", and "optri". 

Despite this PRNG getting its name from the <SearchSpaceIterator> used to generate LCGs, LCGs 
are not required for this PRNG to output values. However, not using a significantly stochastic 
quantity of LCGs, produced from the <SearchSpaceIterator> pass, would mean the <SSINetOp> 
relies much more on its PRNG parameter `prg` to produce reference values for `mdr` and `optri` 
<SSINetNode__TypeLCGNet> types. This basically translates to greater predictability of the 
<SSINetOp> PRNG output, with knowledge of the `prg` output. 

There are two classes of relations connecting the sequence `struct_list` of <SSINetNode__TypeLCGNet>s 
together: 
    - head node (tree identifier) -> directed node tree [`h2tree_map`]
    - tree identifier -> {tree identifiers}. 

These relations determine which <SSINetNode__TypeLCGNet>s get updated, through the course of source 
trees being selected by <SSINetOp> to output values. 

        *Basic Rundown* 
<SSINetOp> starts by setting activation and termination sizes for every non-LCG node. 

Every time <SSINetOp.__next__> is called, PRNG outputs the 0'th value from its primary queue, 
`mainstream_queue`. 

If `mainstream_queue` is empty, <SSINetOp> adds more numbers to queue by this procedure: 
- <SSINetOp> clears its `tmp_queue`. 
- It chooses one of its trees T. 
- It iterates through each node n of T, and adds n's output to `tmp_queue`. 
- For each of the neighboring trees T' connected to T, <SSINetOp> updates the `tmp_queue` 
    of every node of `T'. See the section `*Instantiation Variables*` for the class 
    variables involved in this update process. Also see method<<SSINetOp.distribute_to_connected_trees> 
    and method<SSINetOp.distribute_seq_to_tree>. 
- Additionally, each of these nodes may update their base `struct` by first passing values from their 
  `tmp_queue`s to `base_queue`s, depending on if their counters pass their `activation_size` and 
  `termination_length` parameters. 

    NOTE: there are two primary types for updating: `rapid` and `slow`. See 
    method<SSINetOp.update_node__type_*> for details on this update. 
    This is specified by the class instantiation variable<rapid_update>. 

    The `rapid` update type results in slower number generation. 

""" 
class SSINetOp: 

    def __init__(self,struct_list,h2tree_map,prg,\
        lcg_input_only=0,uniform_io_dist=1,shuffle_dist=0,\
        io_prg=None,rapid_update=False,verbose:bool=False):

        assert len(struct_list) > 3 
        assert len(h2tree_map) > 0
        assert type(prg) in {MethodType,FunctionType}
        assert lcg_input_only in {0,1}
        assert uniform_io_dist in {0,1}
        assert shuffle_dist in {0,1}
        assert type(io_prg) in {MethodType,FunctionType,type(None)}
        assert type(rapid_update) == bool 

        # each element is <SSINetNode__TypeLCGNet> 
        self.struct_list = struct_list
        # `struct_list` index -> Directed Tree 
        # Directed Tree: source node        -> {target nodeset} 
        #               `struct_list` index -> {`struct_list` indices} 
        self.h2tree_map = h2tree_map
        # key of `h2tree_map` -> {other keys of `h2tree_map`} 
        self.edges = dict()
        self.prg = prg 
        self.lcg_input_only = lcg_input_only 
        self.uniform_io_dist = uniform_io_dist 
        self.shuffle_dist = shuffle_dist
        self.io_prg = io_prg
        self.rapid_update = rapid_update 
        self.verbose = verbose 
        
        # storage of values from some source <SSINetNode__TypeLCGNet>; 
        # values can be used for the structures 
        #       `mdr`,`mdrv2`,`optri`. 
        self.mainstream_queue = []
        self.prev_output = None 

        self.tmp_queue = [] 
        self.lcg_queue = []
        self.mo_queue = []

        self.tree_idns = sorted(self.h2tree_map.keys())
        self.t_index = 0

        self.prev = None  

        self.preprocess()
        self.retry_counter = 3 
        return 

    ##-------------------------------------------------------------------

    def preprocess(self): 
        self.set_tree_conn()
        self.set_actterm_sizes()

    def set_tree_conn(self):
        ks = sorted(self.h2tree_map.keys())
        cx = combinations(ks,2)
        q = prg__single_to_int(self.prg)

        for k in ks: 
            self.edges[k] = set() 

        for c in cx: 
            d = int(round(q())) % 2 
            if d:
                if c[0] not in self.edges: 
                    self.edges[c[0]] = set()
                    
                if c[1] not in self.edges:
                    self.edges[c[1]] = set()

                self.edges[c[0]]|= {c[1]}
                self.edges[c[1]]|= {c[0]}

    """
    sets activation and termination sizes for every <SSINetNode__TypeLCGNet> that 
    is not of structure type `lcg`. 
    """
    def set_actterm_sizes(self): 
        for q in self.struct_list: 
            if q.sidn == "lcg": 
                continue 
            q.activation_size = modulo_in_range(int(round(self.prg())),DEFAULT_BASE_QUEUE_ACTIVATION_RANGE) 
            
            qx = [self.prg() for _ in range(q.activation_size)]
            q.update_tmp_queue(qx)

            r = modulo_in_range(self.prg(),\
                DEFAULT_SSINETNODE__TYPE_FITTER_NUM_ITER) 
            q.termination_length = int(ceil(q.activation_size * r)) 
            q.activate_base(self.prg,True)

        return

    ##-------------------------------------------------------------------

    def __next__(self): 

        if self.retry_counter == 0: 
            return 0 

        if len(self.mainstream_queue) > 0: 
            q = self.mainstream_queue.pop(0) 
            self.prev_output = q 
            return q 

        if self.verbose: print("processing trees")

        exclude_trees = set()
        while len(self.mainstream_queue) == 0:
            H = self.choose_tree(exclude_trees)
            if self.verbose: print("- tree {}".format(H))

            # case: no trees produced any values, recursively 
            #       call __next__ 
            if type(H) == type(None): 
                self.retry_counter -= 1 
                return self.__next__() 

            self.process_tree(H) 
            exclude_trees |= {H} 

        q = self.mainstream_queue.pop(0)
        q = self.apply_noise(q) 
        self.prev_output = q 
        self.retry_counter = 3 
        return q 

    """
    50/50 probability of modifying float `q` if `io_prg` is not 
    null. 
    """
    def apply_noise(self,q):
        if type(self.io_prg) == type(None): 
            return q 

        b = int(round(self.io_prg())) % 2 
        if b: 
            return q + self.io_prg() 
        return q 

    def choose_tree(self,exclude_trees=set()):
        TI = [t for t in self.tree_idns if t not in exclude_trees]
        L = len(TI)
        if L == 0: return None 

        i = int(round(self.prg())) % L 
        H = TI[i] 
        return H 

    def process_tree(self,head):
        K = self.tree_to_keys(head,-1)        
        self.clear_tmpqueue()

        # process value for each node 
        for k in K: 
            q = self.process_node(k) 
            if type(q) == type(None): 
                continue 

            if np.isnan(q) or np.isinf(q): 
                continue 

            self.tmp_queue.append(q) 

        if self.verbose: 
            print("\t+ procvec")
            print("\t{}".format(self.lcg_queue))
            print("\t{}".format(self.mo_queue))

        # distribute to other trees 
        self.distribute_to_connected_trees(head)
        self.mainstream_queue.extend(self.tmp_queue)
        return

    def process_node(self,index): 
        q = self.struct_list[index] 

        x = None 
        if q.sidn == "lcg": 
            x = next(q)
            self.lcg_queue.append(x) 
            if self.verbose: print("lcg next: ",x)
            return x 

        if q.sidn in {"mdr","mdrv2"}: 
            x = next(q) 
            if type(x) == type(None): 
                if type(self.prev_output) != type(None): 
                    q.load_first(self.prev_output)
                    x = next(q)
        else:
            x = next(q)

        if self.verbose and q.sidn != "lcg": print("{} next: {}, term={},c={},tmp_q={},base_q={}".\
            format(q.sidn,x,q.termination_length,q.c,len(q.tmp_queue),len(q.base_queue)))

        # case: set new activation and termination length 
        if self.rapid_update: 
            F = self.update_node__type_rapid
        else: 
            F = self.update_node__type_slow

        F(index)         

        if type(x) != type(None): 
            self.mo_queue.append(x)

        return x 

    def update_node__type_slow(self,index): 
        q = self.struct_list[index] 

        if q.c >= q.termination_length:

            q.activation_size = modulo_in_range(int(round(self.prg())),DEFAULT_BASE_QUEUE_ACTIVATION_RANGE) 
            if self.verbose: print("** updating {} @ index={},act={}".format(q.sidn,index,q.activation_size))

            q.c = 0 

            #q.struct = None
            q.base_queue.clear()
            q.base_activated = False
            L = len(q.tmp_queue) 
            if L >= q.activation_size:
                if self.verbose: print("\tbase update")
                q.activate_base(self.prg,self.verbose)

    def update_node__type_rapid(self,index): 
        q = self.struct_list[index] 

        L = len(q.tmp_queue) 
        if L >= q.activation_size:
            if self.verbose: print("\tbase update")
            q.activate_base(self.prg,self.verbose)

        if q.c >= q.termination_length:

            q.activation_size = modulo_in_range(int(round(self.prg())),DEFAULT_BASE_QUEUE_ACTIVATION_RANGE) 
            if self.verbose: print("** updating {} @ index={},act={}".format(q.sidn,index,q.activation_size))

            q.c = 0 

            #q.struct = None
            q.base_queue.clear()
            q.base_activated = False


    #--------- network distribution of output values -----------------------------------

    """
    distributes reference vectors to every neighboring tree connected 
    to `head`. Reference vectors are drawn from either the `lcg_queue` 
    or `tmp_queue`. Each reference vector is appended to its associated 
    node's `tmp_queue`.
    """ 
    def distribute_to_connected_trees(self,head): 
        neighbors = self.edges[head]

        if self.lcg_input_only:
            q = deepcopy(self.lcg_queue)
        else: 
            q = deepcopy(self.tmp_queue)

        prg_ = prg__single_to_int(self.prg) 
        for n in neighbors: 
            if self.shuffle_dist: 
                q = prg_seqsort(q,prg_)
            self.distribute_seq_to_tree(n,q)
        return

    def distribute_seq_to_tree(self,head,L): 

        def D(node): 
            L_ = [self.apply_noise(l) for l in L]

            # case: uniform i/o, send entire sequence to node's tmp queue 
            if self.uniform_io_dist:
                node.update_tmp_queue(deepcopy(L_))
                return 
            
            L2 = [] 
            prg_ = prg__single_to_int(self.prg)
            for l in L_: 
                d = prg_() % 2 
                if d: L2.append(l)
            node.update_tmp_queue(deepcopy(L2)) 

        assert head in self.h2tree_map 
        T = self.h2tree_map[head] 

        if self.verbose: 
            print("distributing to {}:\n\t{}".format(head,L))

        K = self.tree_to_keys(head,1)
        for k in K: 
            node_ = self.struct_list[k] 
            if node_.sidn == "lcg": continue 
            D(node_)
        return 

    def tree_to_keys(self,head,is_bfs): 
        assert head in self.h2tree_map 
        assert is_bfs in {-1,0,1}

        T = self.h2tree_map[head] 
        gd = G2TDecomp(defaultdict(set,T),decomp_rootnodes=[head],child_capacity=1)
        gd.decompose()

        if is_bfs == -1: 
            is_bfs = int(round(self.prg())) % 2 
        tn = gd.decompositions[0]
        
        # collate keys into ordering 
        K = TNode.collate_keys(tn,is_bfs=is_bfs,prg=prg__single_to_int(self.prg))
        return K 

    def clear_tmpqueue(self): 
        self.tmp_queue.clear()
        self.lcg_queue.clear()
        self.mo_queue.clear()
        return 

    #-------------- instantiation ------------------------------------------------ 

    def node_to_activation_count_map(self): 
        return {i:n.uc for (i,n) in enumerate(self.struct_list)} 

    @staticmethod 
    def one_instance(num_nodes,prg,prg2,lcg_input_only=0,uniform_io_dist=1,shuffle_dist=0,\
        prg_io_noise=0,rapid_update=True): 
        assert num_nodes >= 5 
        num_sets = 3 
        var = 0.25 
        num_rounds = 10 
        prt = prg_partition_for_sz__n_rounds(num_nodes,\
            num_sets,prg__single_to_int(prg),var,num_rounds)
        slist = ["lcg","mdr","optri"] 
        
        sls = []

        # make each of the types 
        for (s,p) in zip(slist,prt): 
            sl = SSINetOp.one_instance__slist(s,p,prg)
            sls.extend(sl)

        lx = len(sls)
        q = max([round(lx/3),5])
        q = int(q)
        ssb2n = SSIBatch2Net(sls,q,prg) 
        ssb2n.make_net()

        q1 = sls 
        q2 = ssb2n.h2tree_map
        io_prg = prg if prg_io_noise else None
        return SSINetOp(q1,q2,prg2,lcg_input_only=lcg_input_only,\
            uniform_io_dist=uniform_io_dist,shuffle_dist=shuffle_dist,\
            io_prg=io_prg,rapid_update=rapid_update)

    @staticmethod
    def one_instance__slist(sidn,batch_size,prg):
        # make an LCG from `prg` output 
        ix = [prg() for _ in range(4)] 
        ix = [modulo_in_range(ix_,[1.,101.]) for ix_ in ix]
        prg2 = prg__LCG(ix[0],ix[1],ix[2],ix[3])        
        if sidn == "lcg":
            prg3 = prg__single_to_bounds_outputter(prg,4) 
            param_bounds = prg3()
        elif sidn == "mdr": 
            param_bounds = [0,1]
        else:
            param_bounds = None 

        ssibl = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,batch_size,prg)
        ssibl.instantiate_slist()
        return ssibl.slist 
