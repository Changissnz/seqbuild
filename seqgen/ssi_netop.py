from .ssi_load import * 
from morebs2.numerical_generator import prg_partition_for_sz__n_rounds
from morebs2.g2tdecomp import * 

class SSINetOp:

    def __init__(self,struct_list,h2tree_map,prg,lcg_input_only=-1,uniform_io_dist=1,shuffle_dist=0):
        assert len(struct_list) > 3 
        assert len(h2tree_map) > 0
        assert type(prg) in {MethodType,FunctionType}
        assert lcg_input_only in {-1,0,1}
        assert uniform_io_dist in {0,1}
        assert shuffle_dist in {0,1}

        self.struct_list = struct_list
        self.h2tree_map = h2tree_map
        self.edges = dict()
        self.prg = prg 
        self.lcg_input_only = lcg_input_only 
        self.uniform_io_dist = uniform_io_dist 
        self.shuffle_dist = shuffle_dist

        # storage of values from some source; values can be used for 
        # the structures `mdr`,`mdrv2`,`optri`. 
        self.mainstream_queue = []
        
        self.tmp_queue = [] 
        self.lcg_queue = []
        self.mo_queue = []


        self.tree_idns = sorted(self.h2tree_map.keys())
        self.t_index = 0

        self.prev = None  
        return 

    ##-------------------------------------------------------------------

    def set_tree_conn(self):
        ks = sorted(self.h2tree_map.items())
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

    def set_actsizes(self): 
        for q in self.struct_list: 
            if q.sidn == "lcg": 
                continue 
            q.activation_size = modulo_in_range(int(round(self.prg())),DEFAULT_BASE_QUEUE_ACTIVATION_RANGE) 
        return

    ##-------------------------------------------------------------------

    def __next__(self): 
        return 

    def distribute_to_connected_trees(self,head): 
        T = self.h2tree_map[head] 
        neighbors = self.edges[head]

        if self.lcg_input_only:
            q = self.lcg_queue

        return

    def process_tree(self,head):
        assert head in self.h2tree_map 
        T = self.h2tree_map[head] 

        gd = G2TDecomp(defaultdict(T,set),decomp_nodes=[head],child_capacity=2)
        gd.decompose()

        is_bfs = int(round(self.prg())) % 2 
        tn = gd.decompositions[0]
        
        # collate keys into ordering 
        K = TNode.collate_keys(tn,is_bfs=is_bfs,prg=prg__single_to_int(self.prg))
        self.clear_tmpqueue()

        # process value for each node 
        for k in K: 
            q = self.process_node(k) 

        # distribute to other trees 
        return -1 

    def process_node(self,index): 
        q = self.struct_list[index] 

        x = None 
        if q.sidn == "lcg": 
            x = next(q)
            self.lcg_queue.append(x) 
            return x 

        if q.sidn in {"mdr","mdrv2"}: 
            x = next(q) 
        else:
            x = next(q)

        # case: set new activation and termination length 
        if q.c >= q.termination_length:
            q.activation_size = modulo_in_range(int(round(self.prg())),DEFAULT_BASE_QUEUE_ACTIVATION_RANGE) 
            q.c = 0 

            q.struct = None
            q.base_queue.clear()
            L = len(q.tmp_queue) 
            if L >= q.activation_size:
                q.activate_base(self.prg)
            else: 
                q.base_activated = False

        if type(x) != type(None): 
            self.mo_queue.append(x)

        return x 

    def clear_tmpqueue(self): 
        self.tmp_queue.clear()
        self.lcg_queue.clear()
        self.mo_queue.clear()
        return 

    @staticmethod 
    def one_instance(num_nodes,prg,prg2): 
        assert num_nodes >= 5 

        num_sets = 3 
        var = 0.25 
        num_rounds = 10 
        prt = prg_partition_for_sz__n_rounds(num_nodes,\
            num_sets,prg,var,num_rounds)
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
        return SSINetOp(q1,q2,prg2)

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
