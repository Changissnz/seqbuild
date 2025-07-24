from .n2m import * 
from .vs_fit_op import * 
from .puc import * 
from .ag_ext import * 
from .iseq import * 
from intigers.mod_prng import * 
from morebs2.numerical_generator import prg__constant
from intigers.extraneous import prg__single_to_trinary_vector

DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE = [10,26]
DEFAULT_DX_MULTIPLE_RANGE = [1.0,4]

"""
NOTE: elements of error_vector exist in continuum 
    of negative through positive real numbers. 

    Negative errors are prominently used in difference 
    measures between calculated (expected) and actual. 
"""
def error_improvement_value(error_value):
    q = error_value * -1
    if is_number(q): 
        return round_trinary(q,is_distance_roundtype=False)
    return round_to_trinary_vector(q,is_distance_roundtype=False)

"""
Extension of <IOFit>. Has navigational features for handling 
parameter search space. These features allow this structure to 
improve its solution to the question of the function F, such 
that F(X:input) |--> Y:output. 

Programmed specifically for use with <AffineDelta>. However, other 
class structures can be used as long as they have these functionally 
equivalent methods `fit`,`vectorize`,`update_v2` to that of <AffineDelta>. 

NOTE: to devise alternative fitting structure to the line-based one 
      provided by <AffineDelta>, see those methods in file<vs_fit>. 
"""
class VSSearch(IOFit):

    def __init__(self,x,y,unknown_func,hypdiff_func,madiff_func,\
        prg,depth_rank=10,depth_risk:float=1.0,sol_maxsize:int=1000,\
        is_bfs_queue:bool=False,soln_log_size:int=100): 
        super().__init__(x,y,unknown_func,hypdiff_func,madiff_func)
        assert type(prg) in {FunctionType,MethodType}
        assert type(depth_rank) in {int,np.int32,np.int64} 
        assert depth_rank > 0 
        assert depth_risk >= 0.0 and depth_risk <= 1.0 
        assert type(sol_maxsize) in {int,np.int32,np.int64} 
        assert sol_maxsize > 0 

        self.prg = prg 
        self.depth_rank = depth_rank
        self.depth_risk = depth_risk
        self.sol_maxsize = sol_maxsize
        self.is_bfs_queue = is_bfs_queue
        self.soln_log_size = soln_log_size 
        self.n2mac = None

        self.soln = [] # 
        self.soln_log = [] 
        self.search_queue = []
        self.init_n2m_ac() 

    def init_n2m_ac(self):
        q0,q1 = self.type()

        assert len(q0) == 1
        assert len(q1) == 1

        k0 = q0.pop() 
        k1 = q1.pop() 
        self.n2mac = N2MAutocorrelator((k0,k1),False)
        return
    
    #------------------------- preprocessing methods 
    #------------------------- for initial hypotheses 

    """
    Loads a <HypMach> instance for initial hypotheses on 
    the fitting function F for F(X) = Y. Hypotheses consist 
    of <AffineDelta> instances, calculated via a super-partitioning 
    of the (x_i,y_i) pairs. Super-partitioning uses the 
    condition of being fitted by a common affine function 
    for each of the super-partition sets. 
    """ 
    def preproc(self,ma_order=0): 
        # calculating initial affine hypotheses 
        l = len(self.x) 

        if l < DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE[1]: 
            self.init_HypMach_v2(ma_order,None,\
                None)
            return 
        
        q = modulo_in_range(self.prg(),\
            DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE)
        q = min([q,l // 2])
        assert q > 0 

        index0_vec = []
        index1_vec = []
        indices = [i for i in range(l)] 
        
        index0_vec = prg_choose_n(indices,q,self.prg,is_unique_picker=True) 
        index1_vec = prg_choose_n(indices,q,self.prg,is_unique_picker=True) 

        cx = prg__constant(index1_vec)  
        def index1_function(x): 
            return cx() 

        self.init_HypMach_v2(ma_order,index0_vec,\
            index1_function)
        return 
    
    """
    sets the M+A auxiliary variables of a <HypMach> 
    instance with partially constant values. The constant 
    values are the contribution vector at (0.5,0.5) and 
    a distance vector that consists of values v_i that 
    fall in the range of : (y_i - x_i) x [1.0,4.0]. 

    By this preprocessing scheme, every (x_i,y_i) pair 
    is associated with one affine function. 
    """
    def preproc_v2(self,ma_dim,ma_order): 
        cv = (0.5,0.5) 
        dx = self.default_dx_() 
        self.load_mahyp_auxvar(dx,cv,ma_dim,ma_order) 
        self.init_HypMach() 
    
    """
    generates a distance vector to correspond to 
    the data (`x`,`y`). Distance vector that consists of 
    values v_i that fall in the range of : 
            (y_i - x_i) x [1.0,4.0].
    """
    def default_dx_(self):
        l = len(self.x)
        dx_vec = [] 
        for i in range(l): 
            x,y = self.io_sample(i)
            dx = y - x 

            if is_number(dx):
                r = [int(DEFAULT_DX_MULTIPLE_RANGE[0] * dx),\
                     int(DEFAULT_DX_MULTIPLE_RANGE[1] * dx) + 1] 
                x = modulo_in_range(self.prg(),r) 
                dx_vec.append(x) 
                continue

            r0 = np.array(DEFAULT_DX_MULTIPLE_RANGE[0] * dx,\
                dtype=int) 
            r1 = np.array(DEFAULT_DX_MULTIPLE_RANGE[1] * dx + 1,\
                dtype= int)
            r = np.array([r0,r1]).T 
            dx = [modulo_in_range(self.prg(),r_) for r_ in r] 
            dx_vec.append(dx) 
        return dx_vec 
    
    """
    sets initial hypotheses using <AffineDelta> instances, then 
    clears the variable `mahm`, used to calculate those <AffineDelta>
    instances. 
    """
    def initial_hypotheses(self): 
        q = self.mahm.mhm.info 

        self.hmem = HypMem([],[],mem_type="GHYP") 
        for (i,q_) in enumerate(q):
            h = q_.fit 
            hm0 = self.error_by_hyp(q_.fit) 
            error_term = hm0.c_error(2) 
            gh = GHyp(h,q_,error_term) 
            self.hmem.add(i,gh) 
        self.search_queue = self.hmem.info 
        self.mahm = None 

    #--------------------------------- single hypothesis update methods     

    # TODO: test this 
    """
    applies changes to 1 hypothesis function to produce 
    a sequence of candidate hypothesis functions. These 
    changes are in the form of the `uniform crawl`. The 
    objective is to minimize error. 
    """
    def move_one_hyp__uc(self,unit=10**-1,err_type:int=2,\
        num_attempts:int=1000): 
        assert err_type in {1,2}
        q = self.search_queue.pop(0) 
        
        self.save_soln_to_log()

        # case: depth exceeded
        if q.num_updates >= self.depth_rank: 
            return 

        hm0 = self.error_by_hyp(q.h) 
        vq,p = q.vector_form()
        self.update_soln_set(q,hm0.c_error(2))

        z = np.zeros((len(vq),),dtype=float)

        uc = PUCrawler(z,unit,[len(vq)])
        def ucfunc():
            return next(uc)[0]

        self.move_one__loopty_doo(q,hm0,vq,\
            ucfunc,err_type,num_attempts) 
        return 
    
    """
    moves by auto-correlation 
    """
    def move_one_hyp__ac(self,unit=10**-1,err_type:int=1,\
        num_attempts:int=1000): 
        q = self.search_queue.pop(0) 

        self.save_soln_to_log()
        vq,_ = q.vector_form()

        hm0 = self.error_by_hyp(q.h)
        hm0.condense_error_term(cfunc1=default_cfunc2,\
            cfunc2=default_cfunc1)
        error_value = hm0.c_error(err_type)

        target = error_improvement_value(error_value)
        rdelta = self.rank_xdelta_by_target(target) 
        l = int(round(self.depth_risk * len(rdelta)))
        rdelta = rdelta[:l] 
        rdelta = [r[0] * unit for r in rdelta] 
        if len(rdelta) == 0: 
            return 

        pi = prg__iterable(rdelta,False)
        self.move_one__loopty_doo(q,hm0,vq,pi,err_type,num_attempts)

    """
    moves using output from the pseudo-random number generator 
    `prg`. 
    """
    def move_one_hyp__prg_guided(self,unit=10**-1,err_type:int=1,\
        num_attempts:int=1000): 
        q = self.search_queue.pop(0) 
        self.save_soln_to_log()
        vq,_ = q.vector_form()
        hm0 = self.error_by_hyp(q.h)

        x = prg__single_to_trinary_vector(self.prg,len(vq)) 

        def prg(): 
            return x() * unit 

        self.move_one__loopty_doo(q,hm0,vq,prg,err_type,num_attempts)
        return

    #---------------------------------- methods to aid in adjusting 1 
    #---------------------------------- hypothesis in parameter search. 

    """
    loop process used by method<move_one_hyp__ac>,
    method<move_one_hyp__uc>,method<move_one_hyp__prg_guided>. 
    """
    def move_one__loopty_doo(self,q,hm0,vq,itrtr,err_type,\
        num_attempts:int=1000):
        stat = True 
        c = 0 
        while stat: 
            n = itrtr() 
            c += 1 
            stat = not type(n) == type(None)
            stat = stat and c <= num_attempts
            if not stat: continue             
            q_ = q.update(n,make_copy=True)
            xr,hs = self.cmp_move(hm0,q_.h) 

            xr_ = xr[0] if err_type == 1 else xr[1] 

            # case: there is improvement in cost 
            if xr[1] == 1: 
                q_.num_updates = 0 

            vq2,_ = q_.vector_form()
            stat2 = self.update_soln_set(q_,hs[1].c_error(2))
            if stat2: 
                self.add_back_to_queue(q_) 

            self.n2mac.add_v2(vq,vq2,xr_,True)   

        stat2 = self.update_soln_set(q,hm0.c_error(2))
        if stat2: 
            self.add_back_to_queue(q) 

    def rank_xdelta_by_target(self,target):
        rx = []
        # choose a delta x 
        c = 0 
        for k in self.n2mac.ftable.keys(): 
            kvec = np.array(string_to_vector(k))

            z = None 
            if is_vector(kvec): 
                z = np.zeros((len(kvec),)) 
            else:
                z = 0 
            c += 1 

            j_ = self.n2mac.induce_derivative_v2(\
                z,kvec) 

            if is_number(target) and not is_number(j_):
                j_ = default_cfunc1(j_)

            if is_number(j_) and is_number(target): 
                tv = trinary_diff(j_,target,invertible_weight=2.0) 
            else: 
                tv = trinary_vector_invertible_difference(j_,target,\
                    invertible_weight=2.0) 
                tv = sum(tv)
            rx.append((kvec,tv))
        rx = prg_seqsort_ties(rx,self.prg,vf=lambda x:x[1])
        return rx 
    
    def add_back_to_queue(self,q_): 

        # case: positive feedback, resort to default ordering 
        if q_.num_updates == 0: 
            i = 0 if not self.is_bfs_queue else len(self.search_queue) 
            self.search_queue.insert(i,q_)
            return 
        
        # case: draw number, using prg, that determines  
        #       whether to add `q_` to beginning or end of queue. 
        n = (self.prg() % 1000) / 1000.0
        i = 0 if n < self.depth_risk \
            else len(self.search_queue)
        self.search_queue.insert(i,q_)
        return 

    def cmp_move(self,hm0,h1):
        hm1 = self.error_by_hyp(h1) 
        return hm0.cmp_error(hm1),[hm0,hm1]

    def update_soln_set(self,h,error): 
        # already considered 
        if h.mark: 
            return False 
        
        h.error_term = error 
        h.mark = True 
        if len(self.soln) == 0: 
            self.soln.append((h,error))
            self.soln_size_check()
            return True

        stat = self.already_exists(h)
        if stat: 
            return True  
        
        i = len(self.soln)
        for (j,s) in enumerate(self.soln):
            if s[1] >= error:
                i = j 
                break 

        self.soln.insert(i,(h,error))
        stat = i < self.sol_maxsize 
        self.soln_size_check()
        return stat 

    def already_exists(self,h):
        for s in self.soln: 

            if s[0].hstruct == h.hstruct: 
                return True 
        return False 
    
    def soln_size_check(self):
        while len(self.soln) > self.sol_maxsize: 
            self.soln.pop(-1) 
    

    #--------------------------------------------------
    #---------------------------- methods to measure changes to solution over 
    #---------------------------- the course of search 

    def soln_error_matrix(self):
        q = []
        for x in self.soln_log:
            x2 = [x_[1] for x_ in x]
            q.append(x2) 
        return np.array(q) 

    # TODO: test this. 
    """
    calculates measures on solution log. Outputs 
    4 (4 x 6) matrices. Each of these matrices is 
    of the form: 
    row (0:MIN),(1:MIN),(2:MEAN),(3:VAR) 
    column (0:coverage),(1:uwpd),(2:categorical entropy),(3:min),(4:max),(5:mean). 

    The first two matrices measure the parameter vectors of the 
    solutions, column-wise and row-wise respectively. The last 
    two matrices are likewise measurement values for the error 
    vectors.  
    """
    def measure_soln_log(self):
        if len(self.soln_log) == 0: return None 
        
        def element_to_vector(x):
            x_= []
            x_.extend(x[0])
            x_.append(x[1])
            x_.extend(x[2])
            return x_ 

        def d2matrix(d):
            q = list(d.values())
            q_ = [] 
            for q1 in q: 
                q_.append(element_to_vector(q1)) 
            return np.array(q_) 
        
        def update_accum(acc,s): 
            for i in range(len(acc)): 
                acc[i] += s[i] 

        # measure parameters
        s0accum,s1accum = np.zeros((4,6),dtype=float), \
            np.zeros((4,6),dtype=float)

        for s in self.soln_log:
            q = VSSearch.measure_soln_vector(s,"parameter")
            
            d0 = d2matrix(q[0]) 
            s0 = summarize_matrix(d0) 
            update_accum(s0accum,s0)

            d1 = d2matrix(q[1])
            s1 = summarize_matrix(d1)  
            update_accum(s1accum,s1)

        s0accum /= len(self.soln_log)
        s1accum /= len(self.soln_log) 

        # measure error 
        sm = self.soln_error_matrix()
        d0,d1 = sm.shape 
        prg = prg__iterable(list(sm.flatten()),False)  
        ag = APRNGGaugeV2(prg,(0.,1.),0.5)
        dx2 = ag.measure_matrix_chunk(None,d0,d1,{0,1})

        dxm = d2matrix(dx2[0])         
        s20 = summarize_matrix(dxm) 
        dxm = d2matrix(dx2[1])
        s21 = summarize_matrix(dxm) 

        return (s0accum,s1accum),(s20,s21)

    def save_soln_to_log(self): 
        if len(self.soln) == 0: 
            return 

        q = deepcopy(self.soln)
        self.soln_log.append(q) 

        sz = len(self.soln_log) - self.soln_log_size
        while sz > 0:
            self.soln_log.pop(0)
            sz -= 1


    """
    outputs a measurement for the solution vector `soln`. 

    If `varname` is `error`, outputs 
        (coverage,unidirectional weighted point distance).
    If `varname` is `parameter`, outputs dict 
        D: axis -> index -> 
            [0] (coverage,unidirectional weighted point distance)
            [1] entropy::float 
            [2] (min,max,mean) pairwise difference.
    """
    @staticmethod 
    def measure_soln_vector(soln,varname="error"):
        assert varname in {"error","parameter"} 

        # make the iterator of the terms
        v = []
        d0,d1 = None,None 
        if varname == "error":
            it = [soln_[1] for soln_ in soln]
            prg = prg__iterable(it) 
            ag = APRNGGaugeV2(prg,(0.,1.),0.5)

            return ag.measure_cycle(len(it),\
                term_func=lambda l,l2: type(l) != type(None),\
                auto_frange=True,auto_prange=True,\
                do_cycle_update=True)
        else: 
            itx = [soln_[0].vector_form()[0] for soln_ in soln] 
            it = []
            for itx_ in itx: 
                d0 = len(itx_)
                if type(d1) != type(None): 
                    assert d0 == d1
                d1 = d0 
                it.extend(itx_)
            d1 = len(itx)

            prg = prg__iterable(it) 
            ag = APRNGGaugeV2(prg,(0.,1.),0.5)
            return ag.measure_matrix_chunk(None,d1,d0,{0,1})