from .n2m import * 
from .affine_io_fit import * 
from intigers.tvec import * 
from math import ceil 
from morebs2.numerical_generator import prg__single_to_int,prg_seqsort_ties,prg_decimal

DEFAULT_AFFINE_FIT_SEARCH_GUESS_INCH = 0.1 
DEFAULT_AFFINE_FIT_SEARCH_GUESS_FOOT = 4 + 95 / 100  

DEFAULT_AFFINE_FIT_SPG = 3 
DEFAULT_AFFINE_FIT_CPS = 5 
DEFAULT_AFFINE_FIT_LPG = 12 

# when multiplied with `lifespan_per_guess` for f, minimum frequency for 
# a vector-sign key to be deleted from autocorrelator 
DEFAULT_AFFINE_FIT_KEY_REFRESH_RATIO = 5.0 

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

class AffineDeltaHyp: 

    def __init__(self,ad:AffineDelta,error_term,lifespan,idn): 
        assert type(ad) == AffineDelta
        self.ad = ad 
        self.error_term = error_term
        self.lifespan = lifespan
        self.idn = idn 

        self.is_improvement = False 
        return

    def mark_improvement(self,b:bool): 
        assert type(b) == bool 
        self.is_improvement = b 

    def apply_delta_vector(self,delta_vector,unit): 
        delta = delta_vector * unit
        print("DELTA VECTOR: ",delta_vector)

        i = 0 
        m,a = None,None 

        # get multiple 
        lm = vs_dim(self.ad.m)
        if lm == 0: 
            m = delta[i] 
            i += 1 
        else:
            m = delta[i:i+lm]
            i += lm 
        
        # get additive 
        la = vs_dim(self.ad.a) 
        if la == 0: 
            a = delta[i] 
            i += 1 
        else:
            a = delta[i:i+la]
            i += lm 

        ad = AffineDelta(m,a,self.ad.ma_order)
        return self.ad + ad 

    def __eq__(self,adh): 
        return self.ad == adh.ad 

    def vectorize(self): 
        return self.ad.vectorize() 

    """
    condensed (mean) error 
    """
    def c_error(self): 
        return default_cfunc2(self.error_term)

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
class AffineFitSearch__TypeN2NA(IOAffineFit):

    def __init__(self,x,y,prg,k,ma_dim,ma_order,initial_cv,leap_risk=0.2,\
        is_bfs_queue:bool=False,inch=DEFAULT_AFFINE_FIT_SEARCH_GUESS_INCH,\
        foot=DEFAULT_AFFINE_FIT_SEARCH_GUESS_FOOT,spread_per_guess = DEFAULT_AFFINE_FIT_SPG,\
        candidates_per_spread=DEFAULT_AFFINE_FIT_CPS,lifespan_per_guess = DEFAULT_AFFINE_FIT_LPG,\
        score_improvement_type = "greedy",score_improvement_reference = "local"):

        super().__init__(x,y,prg)

        assert type(k) == int and k > 0 
        assert foot > inch > 0 

        assert 0. <= leap_risk <= 1.0 
        assert score_improvement_type in {"greedy","stochastic"} 
        assert score_improvement_reference in {"lone","local","global"} 

        self.k = k 
        self.ma_dim = ma_dim 
        self.ma_order = ma_order
        self.initial_cv = initial_cv
        self.lrisk = leap_risk
        self.is_bfs_queue = is_bfs_queue

        self.var_dim = None 

        # parameters used to improve on guess 
            # sampling improvement (to test delta of improvement)
        self.inch = inch 
            # attempted improvement 
        self.foot = foot 

            # number of derivatives per guess 
        self.spg = spread_per_guess
            # number of candidates per spread (used for greedy mode) 
        self.cps = candidates_per_spread
            # number of non-improving deltas a guess can take in a row before termination 
        self.lpg = lifespan_per_guess

        self.score_imp_type = score_improvement_type
        self.score_imp_ref = score_improvement_reference

        self.n2mac = None

        self.soln = [] # 

        self.local_best = dict() 
        self.queue = []
        self.preprocess()

        self.finstat = False 

    #---------------------------------- preprocessing 

    def preprocess(self): 
        self.init_hypotheses(self.ma_dim,self.ma_order,self.initial_cv) 
        self.init_n2m_ac()

        self.init_centroids() 

    def var_dim_(self): 
        self.var_dim = np.sum(self.ma_dim)
        if self.ma_dim[0] == 0: self.var_dim += 1 
        if self.ma_dim[1] == 0: self.var_dim += 1 

    def init_n2m_ac(self):
        self.var_dim_() 

        q = self.output_dim 
        if q == 0: q = 1 

        self.n2mac = N2MAutocorrelator((self.var_dim,q),False)
        return

    """
    """
    def init_centroids(self): 

        # calculate five additional centroids 
        extra_ad_guesses = self.averaged_centroid_candidates() 
        
        I = sorted(self.i2f_map.keys()) 
        ads = [self.i2f_map[i] for i in I]
        ads.extend(extra_ad_guesses) 

        error_terms = [(a,self.error_for_sample(a)) for a in ads]

        # rank <AffineDelta>s in ascending order 
        error_terms = sorted(error_terms,key = lambda x: default_cfunc2(x[1])) 
        error_terms = error_terms[:self.k] 

        for (i,et) in enumerate(error_terms): 
            ad,et_ = et 
            H = AffineDeltaHyp(ad,et_,self.lpg,i)
            self.log_soln(H) 
        return 

    """
    used for when user swaps out `Y` for another. 
    """
    def update_soln_err(self): 
        for s in self.soln: 
            E = self.error_for_sample(s.ad) 
            s.error_term = E 
        soln = deepcopy(self.soln) 
        self.soln.clear() 
        self.queue.clear() 
        for s in soln: 
            self.log_soln(s) 


    def averaged_centroid_candidates(self): 
        centroids = [] 
        for _ in range(self.k): 
            c = self.calculate_one_centroid()
            if c not in centroids: 
                centroids.append(c) 
        return centroids 

    def calculate_one_centroid(self): 
        # sample half 
        I = [i for i in range(len(self.X))]
        n = ceil(len(I) / 2) 

        prg_ = prg__single_to_int(self.prg) 
        indices = prg_choose_n(deepcopy(I),n,prg_,True)

        aseq = [self.i2f_map[i] for i in indices]
        return AffineDelta.mean_of_AffineDelta_sequence(aseq,self.ma_dim,self.ma_order) 
    
    #----------------------------------------- main methods 

    def __next__(self): 
        if self.finstat: return 

        if len(self.queue) == 0: 
            self.finstat = True 
            return 

        H = self.queue.pop(0) 
        if self.score_imp_type == "stochastic": 
            self.spread_guess_stochastic(H) 
        else: 
            self.spread_guess__greedy(H)

    def spread_guess_stochastic(self,H): 

        T3 = generate_m_unique_trinary_vectors(self.var_dim,self.spg,self.prg,attempt_ratio=2.0) 

        for T in T3: 
            self.process_one_guess(H,T)

    def spread_guess__greedy(self,H): 

        for _ in range(self.spg): 
            self.spread_guess__greedy_(H) 
        return 

    def spread_guess__greedy_(self,H): 

        T3 = generate_m_unique_trinary_vectors(self.var_dim,self.cps,self.prg,attempt_ratio=2.0) 
        best_tv = self.expected_best_trinary_delta(T3) 
        self.process_one_guess(H,best_tv)

    #----------------------------------------- for guess derivation

    """
    used for greedy approach 
    """
    def expected_best_trinary_delta(self,tvs):  

        tvs_ = [] 
        for t in tvs:
            m = self.n2mac.induce_derivative_v3_(t)
            tvs_.append((t,m))

        # sort with PRNG for tiebreakers 
        vf = lambda x: x[1] 
        Q = prg_seqsort_ties(tvs_,self.prg,vf)
        return Q[0][0] 

    def process_one_guess(self,H,T): 
        # inch 
        H2 = self.cmp_one_guess(H,T,True) 
        self.log_soln(H2) 

        # possible foot 
        d = prg_decimal(self.prg,[0.,1.]) 
        if d > self.lrisk: 
            return  

        H2_ = self.cmp_one_guess(H,T,False) 
        self.log_soln(H2_)

    """
    compares hypothesis H, after application with 
    trinary vector `tvec` to produce H', with 
    reference, specified by `score_imp_ref`. 
    """
    def cmp_one_guess(self,H,tvec,is_inch:bool): 

        # apply the delta 
        X = self.inch if is_inch else self.foot 
        ad_ = H.apply_delta_vector(tvec,X)

        # case: lone 
        if self.score_imp_ref == "lone": 
            err = H.error_term  
        # case: local 
        elif self.score_imp_ref == "local": 
            err = self.best__local_error(H.idn) 
        # case: global 
        else: 
            err = self.best__global_error() 

        err0 = self.error_for_sample(ad_) 

        ##print("E0 {} E1 {}".format(default_cfunc2(err),default_cfunc2(err0)))

        # log into auto-correlator 
        v0 = H.vectorize()[0] 
        v1 = ad_.vectorize()[0] 
        self.log_into_AC(v0,v1,err,err0) 
 
        # compare the two 
        is_improvement = default_cfunc2(err0) < default_cfunc2(err) 

        lpg = self.lpg 
        if not is_improvement: 
            lpg = H.lifespan - 1 

        return AffineDeltaHyp(ad_,err0,lpg,H.idn)

    #--------------------------------------------- retrieving best solutions, based on min. error term 

    def best__global_error(self): 
        assert len(self.soln) > 0
        return self.soln[0].error_term  

    def best__local_error(self,idn): 
        if idn not in self.local_best: 
            return np.ones((self.output_dim,)) * float('inf') 
        return self.local_best[idn].error_term 

    #---------------------------------------------- logging solutions and delta correlations 

    def log_soln(self,H): 

        # case: terminated hypothesis after no improvement 
        if H.lifespan < 0: 
            return 

        # case: non-unique solution 
        if H in self.queue: 
            return 

        # add hypothesis to queue 
        if self.is_bfs_queue:
            self.queue.append(H) 
        else: 
            self.queue.insert(0,H) 

        # rank hypothesis alongside running best global solutions
        # and add it in if qualifies 
        i = len(self.soln) 
        for (j,h_) in enumerate(self.soln): 
            if H.c_error() < h_.c_error(): 
                i = j 
                break 

        self.soln.insert(i,H)

        while len(self.soln) > self.k: 
            self.soln.pop(-1) 

        # check if hypothesis is local best 
        if H.idn not in self.local_best: 
            self.local_best[H.idn] = H 
            return 

        if H.c_error() < self.local_best[H.idn].c_error(): 
            self.local_best[H.idn] = H 

    def error_for_sample(self,s):  
        assert type(s) == AffineDelta

        E = 0 
        if self.output_dim > 0: 
            E = np.zeros((self.output_dim,))

        for (x,y) in zip(self.X,self.Y): 
            y_ = s.fit(x) 
            E += (y - y_)
        return E 

    def log_into_AC(self,x0,x1,e0,e1): 

        # first clear the key at (x1 - x0) if 
        # frequency in `ftable` >= `lpg` * DEFAULT_AFFINE_FIT_KEY_REFRESH_RATIO
        self.n2mac.clear_frequent_key(x0,x1,self.lpg)# * DEFAULT_AFFINE_FIT_KEY_REFRESH_RATIO) 

        self.n2mac.add(x0,x1,e0,e1)