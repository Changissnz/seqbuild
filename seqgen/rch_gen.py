"""
file contains code to extend the <RChainHead> class from the 
library <morebs2>
"""
from intigers.poly_output_fitter_ import * 
from intigers.extraneous import * 
from morebs2.relevance_functions import RCInst,RChainHead
from morebs2.numerical_generator import prg_choose_n


MRIF_VARMAP = {CEPoly:"v",\
    LinCombo:"x"}

DM_FUNC_LIST = [np.dot,mul,safe_div,add,sub]

DEFAULT_RCH_ACCUGEN_RANGE = [-0000,1000]

"""
carries instructions for updating functions <CEPoly>,<LinCombo>, 
set by `base_func`. The `update_freq` is a positive integer that 
specifies how often the class instance accepts a new variable (vector)
to update itself with. 
"""
class MutableRInstFunction:

    def __init__(self,base_func,update_freq:int,idn): 
        assert type(base_func) in MRIF_VARMAP
        assert type(update_freq) == int and update_freq >= 0 
        self.idn = idn 
        self.base_func = base_func 
        self.update_freq = update_freq 
        return

    def dim(self):
        q = getattr(self.base_func,MRIF_VARMAP[type(self.base_func)])
        return len(q) 

    def update(self,new_var): 
        assert type(new_var) == np.ndarray
        if type(self.base_func) == CEPoly: 
            new_var = [(c,i) for (c,i) in enumerate(new_var[::-1])][::-1]
            new_var = np.array(new_var)

        setattr(self.base_func,MRIF_VARMAP[type(self.base_func)],new_var)
        return

    def apply(self,x): 
        y = self.base_func.apply(x)

        if type(y) == np.ndarray: 
            y = y.flatten()
            return safe_npint32_vec(y)
        return safe_npint32_value(y)

"""
<RCHAccuGen> is shorthand for <RChainHead> accumulating generator. 
uses an <RChainHead> `rch` to generate integer values according 
to these specifications: 
- for every `apply` call, the `rch` processes one integer. 
- the `rch` is comprised of a variable number of nodes (type<RCInst>), 
  each with either the <LinCombo.apply> or <CEPoly.apply> function as 
  the variable <RCInst.cf>. The nodes are connected via one uni-directional 
  path. For an `rch` of m nodes, an input value produces (m + 1) values, 
  one at each of the m nodes and the input value itself. 
- for each <RCInst>, the variables `cf` (outer function) and `rf` (value, 
  such as integer, vector) are mutable. 
- the `rf` value (reference value) is used only in the case of <LinCombo>. 
  So for an integer i, the function of the corresponding <RCInst> is
        `cf(dm(i,rf))`; `cf` a <LinCombo.apply> function. 

  The `dm` function is a pairwise function, first argument an integer and the 
  second a vector. See variable `DM_FUNC_LIST`.
- in the processing of integer i, the <RCHAccuGen> accumulates values 
  from `rch.vpath` into its `acc_queue`, the vector of values calculated by 
  each node <RCInst> in the chain. The output from `i` is another integer `j`, 
  the last element from `rch.vpath` at that point in time. 
- to update either `rf` or `cf`, values are drawn from `acc_queue` by use of 
  the `prg` as an index generator. 
- all values stored in `acc_queue` as well as the output value from `apply` 
  are in the bounds of DEFAULT_RCH_ACCUGEN_RANGE. 
"""
class RCHAccuGen: 

    def __init__(self,rch,prg,acc_queue=[],\
        queue_capacity:int=1000,dim_range=(3,8),\
        output_range=DEFAULT_RCH_ACCUGEN_RANGE):

        assert type(acc_queue) == list 
        assert type(queue_capacity) == int and queue_capacity > 1
        self.rch = rch 
        self.prg = prg 
        self.acc_queue = acc_queue 
        self.qcap = queue_capacity
        self.dim_range = dim_range
        self.output_range = DEFAULT_RCH_ACCUGEN_RANGE

        self.mutgen = [set() for _ in range(len(self.rch.s))] 
        self.update_log = defaultdict(defaultdict)
        self.ctr = 0
        return 

    """
    main method
    """
    def apply(self,x): 
        assert is_vector(x) or type(x) in \
            {int,np.int32,np.int64}
        self.rch.apply(x)
        vx = deepcopy(self.rch.vpath)

        for v in self.rch.vpath: 
            if type(v) != np.ndarray: 
                v = safe_npint32_value(v) 
                v = modulo_in_range(v,self.output_range) 
                self.acc_queue.append(v)
            else: 
                v = v.flatten() 
                v = safe_npint32_vec(v) 
                v = [modulo_in_range(v_,self.output_range) \
                    for v_ in v] 
                self.acc_queue.extend(v) 
        self.ctr += 1 
        self.update()
        return self.acc_queue[-1] 

    # TODO: test 
    @staticmethod
    def generate_n_instances(n,prg,num_nodes_range=[1,12],\
        dim_range=[3,6],ufreq_range=[2,9],qcap_range=[100,1001]): 
        
        rgs = [] 
        for i in range(n): 
            num_nodes = modulo_in_range(prg(),num_nodes_range)
            mutrate = modulo_in_range(prg(),[0,1001])
            mutrate = mutrate / 1000.0 
            qcap = modulo_in_range(prg(),qcap_range)
            rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,deepcopy(dim_range),\
                prg,deepcopy(ufreq_range),mutrate,qcap)
            rgs.append(rg) 
        return rgs 

    @staticmethod
    def one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
        ufreq_range,mutrate=0.5,queue_capacity=1000): 

        # declare the RChainHead 
        rch,FX = RCHAccuGen.one_new_RChainHead__v1(\
            num_nodes,dim_range,prg)

        # declare the generator 
        rg = RCHAccuGen(rch,prg,acc_queue=[],\
            queue_capacity=queue_capacity,dim_range=dim_range)

        # choose the candidates to mute
        Q = [i for i in range(num_nodes)]
        n = int(ceil(mutrate * len(Q)))
        Q_ = prg_choose_n(Q,n,prg,is_unique_picker=True)

        for q in Q_: 
            uf = modulo_in_range(prg(),ufreq_range)
            fx = FX[q] 
            mf = MutableRInstFunction(fx,uf,idn=0)
            rg.add_mutable(q,mf)

        return rg 

    @staticmethod
    def one_new_RChainHead__v1(num_nodes,dim_range,prg):
        # declare the RChainHead 
        rch = RChainHead()
        FX = [] 
        for _ in range(num_nodes): 
            rcia,fx = RCHAccuGen.one_new_RCInst_args__v1(dim_range,prg)
            rch.add_node_at(rcia)
            FX.append(fx)
        return rch,FX

    @staticmethod
    def one_new_RCInst_args__v1(dim_range,prg):

        zero = 'r' if prg() % 2 else 'nr' 
        kwargz = [zero]
        FX = None 

        dim = modulo_in_range(prg(),dim_range) 
        # case: referential, uses <LinCombo> 
        if zero == 'r':
            V = safe_npint32__prg_vec(prg,dim)
            kwargz.append(V) 

            fi = prg() % len(DM_FUNC_LIST)
            fx = DM_FUNC_LIST[fi] 
            kwargz.append(fx)

            V2 = safe_npint32__prg_vec(prg,dim + 1)
            lc = LinCombo(V2)

            FX = lc
            kwargz.append(lc.apply)
        # case: non-referential, uses <CEPoly>
        else: 
            V = safe_npint32__prg_vec(prg,dim + 1)
            V2 = [(v,i) for i,v in enumerate(V[::-1])]
            cep = CEPoly(np.array(V2,dtype=np.int32))

            FX = cep
            kwargz.append(cep.apply) 
        return kwargz,FX 

    #------------------- setters/getters for mutable <RChainHead> value 

    def add_mutable(self,rci_index,var): 
        assert rci_index < len(self.mutgen) and \
            rci_index >= 0
        assert type(var) == MutableRInstFunction
        self.mutgen[rci_index] |= {var}


    def tmpset_rch_updatepath(self,rci_index,q,rf_func=lambda x:x): 
        self.rch.s[rci_index].updateInfo = [np.array(q)]
        self.rch.s[rci_index].updateFunc['rf'] = rf_func
        self.rch.s[rci_index].updatePath = {'rf': [0]}
        return

    def delta_on_mutgen(self,rci_index,mg,same_dim:bool=False):
        if len(self.acc_queue) == 0: 
            print("[!] none in queue for op.")
            return 

        assert type(mg) == MutableRInstFunction

        if same_dim:
            d = mg.dim() 
        else:
            d = modulo_in_range(self.prg(),self.dim_range)

        q1 = np.array(prg_choose_n(self.acc_queue,d,self.prg))

        # case: LinCombo, change both 'rf' and 'cf' 
        if type(mg.base_func) == LinCombo: 
            d2 = d - 1
            q2 = np.array(prg_choose_n(self.acc_queue,d2,self.prg))
            self.tmpset_rch_updatepath(rci_index,q2)
            self.rch.s[rci_index].inst_update_var() 

        mg.update(q1) 
        fx = mg.apply
        self.rch.s[rci_index].update_var('cf',fx)
        return

    def mutable2update_list(self): 
        q = [] 
        for (i,x) in enumerate(self.mutgen): 

            for x_ in x: 
                y = self.ctr % x_.update_freq
                if y == 0:
                    q.append((i,x_.idn))

        return q 

    def fetch_mutgen(self,index,idn):
        q = self.mutgen[index]

        for q_ in q: 
            if q_.idn == idn: 
                return q_ 
        return None 

    #--------------------- methods for updating 

    def update(self): 

        # update <RCInst> variables of <RChainHead> 
        ml = self.mutable2update_list()
        for ml_ in ml: 
            self.update_idn(ml_[0],ml_[1])

        # size check for `acc_queue` 
        diff = len(self.acc_queue) - self.qcap
        while diff > 0: 
            self.acc_queue.pop(0) 
            diff -= 1 

        return 

    """
    updates values corresponding to a <MutableRInstFunction>
    """ 
    def update_idn(self,rci_index,var_idn):
        # change appropriate functions/values
        mg = self.fetch_mutgen(rci_index,var_idn)
        self.delta_on_mutgen(rci_index,mg,same_dim=False)

        # edit the update log 
        if rci_index not in self.update_log: 
            self.update_log[rci_index] = defaultdict(int) 
        self.update_log[rci_index][var_idn] += 1 
        return 

    def __next__(self):
        return -1

    def change_dim(self): 
        return -1 