from .ssi_load import * 
from morebs2.numerical_generator import prg_partition_for_sz__n_rounds
from morebs2.g2tdecomp import * 

class SSINetOp:

    def __init__(self,struct_list,h2tree_map,prg,\
        lcg_input_only=0,uniform_io_dist=1,shuffle_dist=0,\
        io_prg=None,verbose:bool=False):

        assert len(struct_list) > 3 
        assert len(h2tree_map) > 0
        assert type(prg) in {MethodType,FunctionType}
        assert lcg_input_only in {0,1}
        assert uniform_io_dist in {0,1}
        assert shuffle_dist in {0,1}
        assert type(io_prg) in {MethodType,FunctionType,type(None)}

        self.struct_list = struct_list
        self.h2tree_map = h2tree_map
        self.edges = dict()
        self.prg = prg 
        self.lcg_input_only = lcg_input_only 
        self.uniform_io_dist = uniform_io_dist 
        self.shuffle_dist = shuffle_dist
        self.io_prg = io_prg

        self.verbose = verbose 
        
        # storage of values from some source; values can be used for 
        # the structures `mdr`,`mdrv2`,`optri`. 
        self.mainstream_queue = []
        self.prev_output = None 

        self.tmp_queue = [] 
        self.lcg_queue = []
        self.mo_queue = []

        self.tree_idns = sorted(self.h2tree_map.keys())
        self.t_index = 0

        self.prev = None  

        self.preprocess()
        return 

    ##-------------------------------------------------------------------

    def preprocess(self): 
        self.set_tree_conn()
        self.set_actterm_sizes()

    def set_tree_conn(self):
        ks = sorted(self.h2tree_map.keys())
        cx = combinations(ks,2)
        q = prg__single_to_int(self.prg)

        for c in cx: 
            d = int(round(q())) % 2 
            if d:
                if c[0] not in self.edges: 
                    self.edges[c[0]] = set()
                    
                if c[1] not in self.edges:
                    self.edges[c[1]] = set()

                self.edges[c[0]]|= {c[1]}
                self.edges[c[1]]|= {c[0]}

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
            q.activate_base(self.prg)

        return

    ##-------------------------------------------------------------------

    def __next__(self): 

        if len(self.mainstream_queue) > 0: 
            q = self.mainstream_queue.pop(0) 
            self.prev_output = q 
            return q 

        if self.verbose: print("processing trees")

        exclude_trees = set()
        while len(self.mainstream_queue) == 0:
            H = self.choose_tree(exclude_trees)
            if self.verbose: print("- tree {}".format(H))
            assert type(H) != type(None)
            
            self.process_tree(H) 
            exclude_trees |= {H} 

        q = self.mainstream_queue.pop(0)
        q = self.apply_noise(q) 
        self.prev_output = q 
        return q 

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
            if type(q) != type(None): 
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
            if self.verbose: print("{} next: {}".format(q.sidn,x))
        else:
            x = next(q)
            if self.verbose: print("{} next: {}".format(q.sidn,x))

        # case: set new activation and termination length 
        if q.c >= q.termination_length:
            if self.verbose: print("** updating {} @ index={}".format(q.sidn,index))

            q.activation_size = modulo_in_range(int(round(self.prg())),DEFAULT_BASE_QUEUE_ACTIVATION_RANGE) 
            q.c = 0 

            #q.struct = None
            q.base_queue.clear()
            q.base_activated = False
            L = len(q.tmp_queue) 
            if L >= q.activation_size:
                if self.verbose: print("\tbase update")
                q.activate_base(self.prg,self.verbose)

        if type(x) != type(None): 
            self.mo_queue.append(x)

        return x 


    #--------- network distribution of output values -----------------------------------

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
            if self.uniform_io_dist:
                node.update_tmp_queue(deepcopy(L_))
                return 
            
            L_ = [] 
            prg_ = prg__single_to_int(self.prg)
            for l in L: 
                d = prg_() % 2 
                if d: L_.append(l)
            node.update_tmp_queue(deepcopy(L_))

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

    @staticmethod 
    def one_instance(num_nodes,prg,prg2,lcg_input_only=0,uniform_io_dist=1,shuffle_dist=0,\
        prg_io_noise=0): 
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
            uniform_io_dist=uniform_io_dist,shuffle_dist=shuffle_dist,io_prg=io_prg)

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
