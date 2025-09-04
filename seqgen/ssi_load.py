from morebs2.search_space_iterator import * 
from intigers.lcg_v3 import *  
from intigers.prng_pw_op import * 
from intigers.extraneous import * 

SSIBL_STRUCTURES = ["lcg","mdr","mdrv2","optri"]

SSIBL_PARAMETER_MAP = {'lcg':[float,float,float,float],\
        'mdr':[list,bool],\
        'optri':[list,"op_f","op_b"]}

DEFAULT_SSINET_HOP_RANGE = [4,10] 

class SSINetNode__TypeLCGNet: 

    def __init__(self,structure_idn,prg,op_pair=None):
        assert structure_idn in SSIBL_STRUCTURES
        assert type(prg) in {MethodType,FunctionType,type(None)}
        self.sidn = structure_idn
        self.prg = prg 

        assert type(op_pair) in {type(None),tuple}
        self.op_pair = op_pair
        self.out_queue = [] 
        return 

class SSIBatchLoader__TypeLCGNet: 

    """
    param_bounds := bounds vector (if structure is lcg) | ordered binary range (i.e. 00,01,10,11)
    """
    def __init__(self,structure_idn,param_bounds,max_batch_size,aux_prg): 
        assert structure_idn in SSIBL_STRUCTURES
        self.sidn = structure_idn
        self.check_param_bounds(param_bounds)

        self.index = 0 # used for ordered binary range  
        assert max_batch_size > 0 
        self.max_batch_size = max_batch_size

        assert type(aux_prg) in {MethodType,FunctionType}
        self.aux_prg = aux_prg 
        # list of structures 
        self.slist = [] 

        if self.sidn == "lcg": 
             self.load_SSI()

        return 

    def check_param_bounds(self,param_bounds):
        if self.sidn == "lcg":
            stat0 = is_bounds_vector(param_bounds)
            assert stat0 
            assert np.all(param_bounds[:,1] >= param_bounds[:,0]) 
            self.param_bounds = param_bounds

        elif "mdr" in self.sidn:
            assert len(param_bounds) == 2 
            assert min(param_bounds) in {0,1}
            assert max(param_bounds) in {0,1}
            assert param_bounds[0] <= param_bounds[1]
            self.param_bounds = param_bounds

        elif self.sidn == "optri": 
            self.param_bounds = None  
        return -1 

    def instantiate_slist(self):
        while self.max_batch_size > 0: 
            s = self.one_struct() 
            if type(s) == type(None): 
                break 
            self.slist.append(s)
            self.max_batch_size -= 1
             
        return 

    def one_struct(self):
        if self.sidn == "lcg": 
            if self.ssi.finished(): 
                return None 
            rs = next(self.ssi) 
            prg = prg__LCG(rs[0],rs[1],rs[2],rs[3]) 
            return SSINetNode__TypeLCGNet(self.sidn,prg,op_pair=None)

        elif self.sidn == "mdr":
            return SSINetNode__TypeLCGNet(self.sidn,None,op_pair=None)

        else:
            o1,o2 = self.one_operator_pair()
            return SSINetNode__TypeLCGNet(self.sidn,None,op_pair=(o1,o2))

    def load_SSI(self):
        assert self.sidn == "lcg" 

        aux_prg_ = prg__single_to_int(self.aux_prg)
        column_order = prg_seqsort([0,1,3,2],aux_prg_)
        ssihops = [int(modulo_in_range(self.aux_prg(),DEFAULT_SSINET_HOP_RANGE)) for _ in range(4)]

        for i,r in enumerate(self.param_bounds): 
            if abs(r[1] -r[0]) < 0.05: 
                ssihops[i] = 2
        ssihops = np.array(ssihops)
        self.ssi = SearchSpaceIterator(self.param_bounds, startPoint=deepcopy(self.param_bounds[:,0]), \
            columnOrder=column_order, SSIHop = ssihops,cycleOn = False, cycleIs = 0)

    def one_operator_pair(self):
        base_ops = DEFAULT_PAIRWISE_OPS + [mod]

        q = prg__one_weighted_pairwise_operator(self.aux_prg,base_ops,\
            base_ops,weight_range=(-2.222,2.555))

        q1 = prg__one_weighted_pairwise_operator(self.aux_prg,base_ops,\
            base_ops,weight_range=(-2.222,2.555))

        return q,q1 


class SSIBatch2Net:

    def __init__(self,structure_instances,num_trees,prg): 
        assert len(structure_instances) > 0 
        assert type(num_trees) == int 
        assert num_trees <= len(structure_instances) 
        assert type(prg) in {MethodType,FunctionType}

        self.si = structure_instances
        self.num_trees = num_trees 
        self.prg = prg 
        self.h2tree_map = dict()
        return 

    def make_net(self): 
        for _ in range(self.num_trees): 
            h,t = self.one_tree() 
            self.h2tree_map[h] = t  
        return 

    """
    """ 
    def one_tree(self): 
        qi_ = set([i for i in range(len(self.si))]) 
        qi = qi_ - set(self.h2tree_map.keys()) 
        assert len(qi) > 0
        qi = sorted(qi)

        i = int(self.prg() % len(qi))

        # choose head 
        head = qi[i] 
        queue = [head]
        cache = []
        qi_ = sorted(qi_ - {head})

        tx = dict()
        while len(qi_) > 0 and len(queue) > 0: 
            hx = queue.pop(0)
            cache.append(hx)

            # choose neighbors 
            nn = int(modulo_in_range(self.prg(),(1,len(qi_) + 1)))
            neighbors = set()
            for _ in range(nn): 
                i = int(self.prg() % len(qi_)) 
                n = qi_.pop(i)
                neighbors = neighbors | {n}
            tx[hx] = neighbors 

            queue.extend(sorted(neighbors))
            qi_ = [q for q in qi_ if (q not in neighbors and q != hx)] 

        return head,tx 
    