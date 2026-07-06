from desi.multi_metric import * 
from desi.seqcov_perm import * 
from collections import deque 
from morebs2.numerical_generator import prg_decimal

DEFAULT_AGV2GG_NGRAM_RANGE = [4,9] 
DEFAULT_AGV2GG_SUBRANGE_SIZE_RANGE = [1,5] 
DEFAULT_AGV2GG_SUBRANGE_RESET_TIME_RANGE = [11,200]
DEFAULT_AGV2GG_CATEGORY_SIZE_RANGE = [3,11] 

DEFAULT_AGV2GG_AUTODIST_COVERAGE_ATTEMPTS = 3 
DEFAULT_AGV2GG_COVERAGE_RADIUS = round(49/101,5) 

DEFAULT_AGV2GG_AUTODIST_UWPD_ATTEMPTS = 8

def uwpd_of_modded_sequence(S,m):
    assert m != 0 
    S_ = np.array([s_ % m for s_ in S]) 
    return uwpd(S_,pairwise_op=lambda x1,x2: np.abs(x2 - x1),\
        accum_op=lambda x1,x2: x1 + x2)

def uwpdcov_of_modded_sequence(S,m): 
    pdseq = uwpd_of_modded_sequence(S,m)
    mseq = np.array([s_ % m for s_ in S]) 
    covseq = coverage_of_sequence(mseq,[0,m],max_radius=0.5)
    return pdseq,covseq 

"""

base_output_span := int, positive value specifying the length of the vector from the 
                `base_prg` to use.
density_measure := int, positive value specifying the maximum length of the category 
                vector used by this PRNG for output to satisfy metric objectives.
ngram_length := int, positive value at most equal to `base_output_span`, specifying 
            the n-gram length (subvectoring) for the current `base_vec`. 
allow_subrange_drift := bool, if set to True, the PRNG resorts to outputting numbers in 
        subranges of the `super_range`. 
"""
class GaugeGuidedGen: 

    def __init__(self,base_prg,\
        base_output_span,density_measure:int,super_range,\
        ngram_length:int=None,target_measure="cov",allow_subrange_drift:bool=False):

        assert type(base_prg) in {MethodType,FunctionType}
        assert base_output_span[0] > 3 and is_valid_range(base_output_span,True,False), \
            "base output span has to be at least 3"

        assert density_measure >= 5 and type(density_measure) == int, "density measure has to be at least 5"
        assert is_valid_range(super_range,True,False) or is_valid_range(super_range,False,False)
        
        if type(ngram_length) != type(None): 
            assert type(ngram_length) == int and ngram_length > 3, "ngram length has to be at least 3" 
        else: 
            ngram_length = modulo_in_range(int(base_prg()),DEFAULT_AGV2GG_NGRAM_RANGE)

        assert target_measure in {"cov","uwpd"}
        assert type(allow_subrange_drift) == bool 

        self.base_prg = base_prg 
        self.base_output_span = base_output_span
        self.density = density_measure
        self.super_range = super_range 
        self.ngram_length = ngram_length
        self.target_measure = target_measure 
        self.allow_subrange_drift = allow_subrange_drift

        self.queue = deque() 
        self.mvec = deque()  
        self.catvec = deque() 

        self.cat_size = None 

        self.prev_vec = None 
        self.base_vec = np.array([])  

        self.active_ranges = [] 

        self.sr_counter = 0 
        self.sr_counter_max = None  
        self.set_range() 

    def set_range(self): 
        self.active_ranges.clear()

        self.sr_counter_max = modulo_in_range(int(self.base_prg()),\
            DEFAULT_AGV2GG_SUBRANGE_RESET_TIME_RANGE)

        if not self.allow_subrange_drift: 
            self.active_ranges = [self.super_range] 
            return 

        n = modulo_in_range(int(self.base_prg()),DEFAULT_AGV2GG_SUBRANGE_SIZE_RANGE)
        prg_ = wrap_ranged_modulo_over_generator(self.base_prg,self.super_range) 
        S = sorted(prg_unique_sequence(prg_,n * 2)) 

        for i in range(n): 
            self.active_ranges.append(S[i*2:(i*2)+2]) 
        return 

    def __next__(self): 

        if len(self.queue) == 0: 
            V = self.next__base()
            self.base_vec = deepcopy(V)  
            self.queue.extend(V) 

        self.maintain_density() 
        self.sr_counter += 1 
        if self.sr_counter >= self.sr_counter_max: 
            self.set_range() 
            self.sr_counter = 0 

        return self.queue.popleft() 

    def maintain_density(self):
        while len(self.mvec) > self.density: 
            self.mvec.popleft() 
        
        while len(self.catvec) > self.density: 
            self.catvec.popleft()

    def next__base(self): 
        l = modulo_in_range(int(self.base_prg()),self.base_output_span)
        x = prg__single_to_nvec(self.base_prg,l)() 

        self.next_base_vec(x)

        self.cat_size = modulo_in_range(int(self.base_prg()),\
            DEFAULT_AGV2GG_CATEGORY_SIZE_RANGE)

        return self.auto_distribute() 

    def next_base_vec(self,x): 

        if type(self.prev_vec) == type(None) or len(self.prev_vec) == 0: 
            self.prev_vec = deepcopy(self.base_vec)
        else: 
            self.prev_vec = modulated_vec_op(self.prev_vec,self.base_vec,add) 
            self.prev_vec = modulo_in_range(self.prev_vec,self.super_range) 


        self.base_vec = modulo_in_range(x,self.super_range)
        self.base_vec = np.round(self.base_vec,5)

        if len(set(self.base_vec)) < 2: 
            s = prg_unique_sequence(self.base_prg,len(self.base_vec))
            self.base_vec += s 


    def auto_distribute(self): 

        # case: none to compare 
        if len(self.mvec) == 0: 
            m0,m1 = self.measure(self.base_vec)
            self.mvec.append((m0,m1))
            
            c = self.classify(m0,m1)
            self.catvec.append(c) 
            return self.base_vec 
        
        # case: get the best scoring derivation 
        V,m0,m1,c = self.best_distribution() 

        self.mvec.append((m0,m1))
        self.catvec.append(c) 
        return V 

    def best_distribution(self):  

        V_ = self.base_vec  

        # add this one to previous vector 
        if type(self.prev_vec) != type(None): 
            V_ = modulated_vec_op(V_,self.base_vec,add) 
            V_ = modulo_in_range(V_,self.super_range)

        if self.target_measure == "cov": 
            return self.best_distribution__coverage(V_) 
        return self.best_distribution__uwpd(V_)        

    def best_distribution__coverage(self,V_): 
        m0_,m1_,c_,s_ = None,None,-1,-float('inf')
        V = None 

        for i in range(1,DEFAULT_AGV2GG_AUTODIST_COVERAGE_ATTEMPTS + 1): 
            r = [(i - 1) / DEFAULT_AGV2GG_AUTODIST_COVERAGE_ATTEMPTS, i / DEFAULT_AGV2GG_AUTODIST_COVERAGE_ATTEMPTS] 
            d = prg_decimal(self.base_prg,r)
            V = self.permute_coverage(V_,d)
            m0,m1,c,s = self.compare_with_previous(V)
            if s > s_: 
                m0_,m1_,c_,s_ = m0,m1,c,s 
                V_ = V 

        return V_,m0_,m1_,c_ 

    def best_distribution__uwpd(self,V_): 

        l = modulo_in_range(int(self.base_prg()),self.base_output_span)
        A = prg__single_to_nvec(self.base_prg,len(self.base_vec))() 

        Q = add_trinary_vector_based_noise(V_,A,self.base_prg,\
            DEFAULT_AGV2GG_AUTODIST_UWPD_ATTEMPTS)

        m0_,m1_,c_,s_ = None,None,-1,-float('inf')
        V = None 
        for q in Q: 
            m0,m1,c,s = self.compare_with_previous(q) 
            if s > s_: 
                m0_,m1_,c_,s_ = m0,m1,c,s 
                V_ = q  

        return V_,m0_,m1_,c_ 

    def permute_coverage(self,V,coverage_delta): 

        sub_range = self.active_ranges[int(self.base_prg()) % len(self.active_ranges)]
        V = modulo_in_range(V,sub_range) 

        scp = SeqCoveragePermuter(V,coverage_delta,\
            DEFAULT_AGV2GG_COVERAGE_RADIUS,sub_range,self.base_prg,delta_type=1)
        scp.set_partition(self.cat_size) 
        return scp.apply() 

    def permute_uwpd(self,V): 
        sub_range = self.active_ranges[int(self.base_prg()) % len(self.active_ranges)]
        sup = SeqUWPDPermuter(V,uwpd_delta,sub_range,self.base_prg)
        return sup.apply()

    def compare_with_previous(self,V): 
        if len(self.mvec) == 0: return None,None,None 

        m0,m1 = self.measure(V) 
        c = self.classify(m0,m1)
        s = index_based_weighted_distance(c,self.catvec)
        return m0,m1,c,s 

    def classify(self,m0,m1): 
        m = m0 
        if self.target_measure == "cov":
            m = m0
        else:
            m = m1 

        return std_classify_one_value(m,1.0/self.cat_size,0.0)

    """
    return:
    - mean of coverage, mean of u.w.p.d. 
    """
    def measure(self,V):  
        mm = MultiMetric(V)
        smry = mm.summarize(self.ngram_length,False)

        q0 = smry[0][:,0]
        q1 = smry[0][:,1] 

        return np.mean(q0),np.mean(q1) 

