from morebs2.search_space_iterator import * 
from morebs2.matrix_methods import is_2dmatrix 

from intigers.lcg_v3 import *  
from intigers.prng_pw_op import * 
from intigers.extraneous import * 

SSIBL_STRUCTURES = ["lcg","mdr","mdrv2","optri"]

SSIBL_PARAMETER_MAP = {'lcg':[float,float,float,float],\
        'mdr':[list,bool],\
        'optri':[list,"op_f","op_b"]}

DEFAULT_SSINET_HOP_RANGE = (4,10) 
DEFAULT_BASE_QUEUE_ACTIVATION_RANGE = (5,24)    
DEFAULT_SSINETNODE__TYPE_FITTER_NUM_ITER = (0.5,4.1)  
DEFAULT_SSINETNODE_TMPCACHE_RATIO = 4.9  

def npint32__mod__Q(v):
    r = (-5000,5000) 
    if v >= r[0] and v <= r[1]: return np.int32(v) 
    v_ = modulo_in_range(v,r) 
    return np.int32(v_)

def index1_to_index2(i,dim2d):
    x = 0 
    i2 = [0,0]
    while x < i:  
        x += dim2d[1] 

        if x <= i: 
            i2[0] += 1  
    
    r = i - i2[0] * dim2d[1]
    i2[1] = r 
    return i2

class OpTriGenLite: 

    def __init__(self,iseq,forward_func,backward_func):
        assert type(iseq) == IntSeq 
        self.ffunc = forward_func
        self.bfunc = backward_func
        self.to_matrix(iseq) 
        self.i = 0 
        return 

    def to_matrix(self,iseq): 
        ot1 = iseq.optri(self.bfunc,float) 
        iseq2 = IntSeq(ot1[0])
        ot2 = iseq2.optri(self.ffunc,float)

        ot2 = np.flip(ot2,1)
        ot2 = np.flip(ot2,0)
        j = 1
        for i in range(1,len(ot1)): 
            x = ot2[i-1][:j] 
            ot1[i,:j] = x 
            j += 1  
        self.m = ot1 

    def __next__(self):
        qi = index1_to_index2(self.i,self.m.shape)
        self.i = (self.i + 1) % (self.m.shape[0] * self.m.shape[1]) 
        return self.m[qi[0],qi[1]]

class SSINetNode__TypeLCGNet: 

    def __init__(self,structure_idn,prg,op_pair=None,exclude_neg=None,\
        tmpcache_ratio= DEFAULT_SSINETNODE_TMPCACHE_RATIO):
        assert structure_idn in SSIBL_STRUCTURES
        assert type(prg) in {MethodType,FunctionType,type(None)}
        self.sidn = structure_idn

        # used for `lcg` 
        self.prg = prg 

        assert type(op_pair) in {type(None),tuple}
        if type(op_pair) == tuple: 
            assert len(op_pair) == 2 
        assert tmpcache_ratio >= 1.0 

        # used for `mdr` 
        self.exclude_neg = exclude_neg
        self.tmpcache_ratio = tmpcache_ratio

        # used for `optri` 
        self.op_pair = op_pair 

        # used for `optri`, `mdr` 
        self.base_queue = [] 
        self.base_activated = False 
        self.tmp_queue = [] 

        self.struct = None 
        # use for `mdr` `
        self.rvalues = [] 

        self.activation_size = None
        self.termination_length = None 

        self.c = 0 
        return 

    def __next__(self):
        if self.sidn == "lcg": 
            return self.prg()
        
        if not self.base_activated: 
            return None 

        if self.sidn in {"mdr","mdrv2"}: 
            if len(self.rvalues) == 0: 
                q = None  
                self.c -= 1 
            else:
                q = self.rvalues.pop(0)
        else:
            q = next(self.struct)
            
        self.c += 1 
        return q 

    """

    """
    def set_actsize(self, activation_size): 
        assert "mdr" in self.sidn or "optri" in self.sidn
        assert type(activation_size) in {int,np.int32} 
        assert DEFAULT_BASE_QUEUE_ACTIVATION_RANGE[1] >= \
            activation_size >= DEFAULT_BASE_QUEUE_ACTIVATION_RANGE[0]
        self.activation_size = activation_size
        return

    def activate_base(self,aux_prg,verbose=0): 
        assert type(aux_prg) in {MethodType,FunctionType}
        
        self.base_activated = True 
        x = self.tmp_queue[:self.activation_size]
        self.base_queue = x 
        self.tmp_queue = self.tmp_queue[self.activation_size:]

        q = self.activate_mdr(aux_prg,verbose)
        if q == False:
            q = self.activate_optri(aux_prg,verbose) 
        self.struct = q 
        return 

    def activate_mdr(self,aux_prg,verbose):    
        if "mdr" not in self.sidn: 
            return False 

        BQ = [npint32__mod__Q(b) for b in self.base_queue]

        mdr = None
        if verbose: 
            print("= activating {} w/ \n{}".format(self.sidn,BQ))
        
        if self.sidn == "mdr":
            m = ModuloDecomp(IntSeq(BQ))
            m.merge(bool(self.exclude_neg)) 
            mdr = ModuloDecompRepr(m,1)
        else: 
            m = ModuloDecomp(IntSeq(BQ),\
                self.exclude_neg)
            mdr = ModuloDecompRepr(m,2) 

        L = len(self.base_queue)
        r = modulo_in_range(aux_prg(),\
            DEFAULT_SSINETNODE__TYPE_FITTER_NUM_ITER) 
        L_ = int(ceil(L * r)) 
        self.termination_length = L_ 
        return mdr

    def activate_optri(self,aux_prg,verbose): 
        L = len(self.base_queue) - 1 
        L = L ** 2 
        r = modulo_in_range(aux_prg(),\
            DEFAULT_SSINETNODE__TYPE_FITTER_NUM_ITER) 
        L_ = int(ceil(L * r)) 
        self.termination_length = L_

        if verbose: 
            print("= activating {} w/ \n{}".format(self.sidn,self.base_queue))

        B = [npint32__mod__Q(b) for b in self.base_queue]
        iseq = IntSeq(B) 
        return OpTriGenLite(iseq,self.op_pair[0],self.op_pair[1])

    def update_tmp_queue(self,L): 
        assert type(L) == list
        assert self.sidn != "lcg"

        self.tmp_queue.extend(L) 

        limit = int(round(self.activation_size * self.tmpcache_ratio)) 
        excess = len(self.tmp_queue) - limit 
        while excess > 0: 
            self.tmp_queue.pop(0)
            excess-= 1
        return

    def load_first(self,v): 
        assert is_number(v) 
        self.struct.reset_first(v,False)
        self.rvalues = self.struct.reconstruct()
        return

class SSIBatchLoader__TypeLCGNet: 

    """
    param_bounds := bounds vector (if structure is lcg) | ordered binary range (i.e. 00,01,10,11)
    """
    def __init__(self,structure_idn,param_bounds,max_batch_size,aux_prg,ssi_hr=DEFAULT_SSINET_HOP_RANGE): 
        assert structure_idn in SSIBL_STRUCTURES
        self.sidn = structure_idn
        self.check_param_bounds(param_bounds)

        self.index = 0 # used for ordered binary range  
        assert max_batch_size > 0 
        self.max_batch_size = max_batch_size

        assert type(aux_prg) in {MethodType,FunctionType}
        self.aux_prg = aux_prg 

        assert is_valid_range(ssi_hr,True,False)
        self.ssi_hr = ssi_hr

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
            b = int(round(modulo_in_range(self.aux_prg(),[0.,1.])))
            return SSINetNode__TypeLCGNet(self.sidn,None,op_pair=None,exclude_neg=bool(b))

        else:
            o1,o2 = self.one_operator_pair()
            return SSINetNode__TypeLCGNet(self.sidn,None,op_pair=(o1,o2))

    def load_SSI(self):
        assert self.sidn == "lcg" 

        aux_prg_ = prg__single_to_int(self.aux_prg)
        column_order = prg_seqsort([0,1,3,2],aux_prg_)
        ssihops = [int(modulo_in_range(self.aux_prg(),self.ssi_hr)) for _ in range(4)]

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
    