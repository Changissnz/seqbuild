from mini_dm.ngram_ext import * 
from morebs2.seq_repr import * 
from types import MethodType,FunctionType
from mini_dm.iseq import stdop_vec

"""
Used for Cycle Test, focused on qualities of forward permutations. Measures are 
frequency-based on n-gram prefixes to m-gram suffixes of the next `num_iter` samples 
from `prg`. 

Important metrics per prefix: 
[0] mcs_k, MCS k-complexity of the ordered concatenation of suffixes. 
[1] cy_diff, mean cyclical kernel difference of ordered concatenation of suffixes. 
[2] I_diff, set size for difference vector of prefix indices of occurrence.
[3] var_norm, normalized variance of suffix frequency. 

For metrics [0] + [1], there are an additional four categories of scores: 
    seq -> ordered concatenation of literal suffixes,
    set -> ordered concatenation of set form of suffixes,
    valfreq -> ordered concatention of (suffix label,suffix frequency) pairs per suffix, 
    l-seq -> ordered concatenation of labelled suffixes. 
"""
class CycleMetric: 

    def __init__(self,prg,num_iter:int,pre_size:int,post_size:int,cast_func=cr):
        assert type(prg) in {MethodType,FunctionType} 
        assert num_iter > 0 and type(num_iter) == int  

        self.prg = prg 
        self.num_iter = num_iter 
        self.pre_size = pre_size 
        self.post_size = post_size 
        self.cast_func = cast_func 

        self.S = [] 
        self.M_ = None  
        self.M = [] 

        #   prefix --> measure type --> (mcs_k, cy_diff)
        self.prefix2mcs_map = dict()

        #   prefix --> (I_diff, var_norm)
        self.prefix2score_map = dict() 

        self.preproc() 

    def preproc(self): 
        self.load_next_batch() 
    
    """
    call this to update reference sequence to another 
    one, before calling method<process>. 
    """
    def load_next_batch(self,num_iter=None):
        if type(num_iter) != type(None): 
            assert type(num_iter) == int and num_iter > 0 
            self.num_iter = num_iter 

        self.S.clear() 
        self.S = [self.cast_func(self.prg()) for _ in range(self.num_iter)]  
        
        if type(self.M_) != type(None): 
            self.M.append(self.M_)

        self.M_ = SequencePrePostMap(self.S,self.pre_size,self.post_size) 
        self.M_.process()

        self.prefix2score_map.clear()
        for k in self.M_.M.keys(): 
            self.prefix2score_map[k] = dict() 
    
    """
    main method
    """
    def process(self): 

        prefixes = set(self.M_.M.keys())
        S = dict() 
        print("num prefixes: ",len(prefixes))
        for p in prefixes: 
            print("processing: ",p)
            P,d0,v0 = self.process_prefix(p) 
            self.prefix2mcs_map[p] = P 
            self.prefix2score_map[p] = (d0,v0) 
        print("")
        return

    """
    [0] dict, measure type --> mean over prefixes
    [1] mean of I_diff over prefixes 
    [2] mean of normed variance over prefixes 
    """
    def summarize(self): 
        D = dict() 
        for m in DEFAULT_CYMETRIC_SUFFIX_MEASURE_TYPES: 
            S = self.fetch_scores_by_type(m) 
            M = np.array(list(S.values())) 
            q = np.mean(M,axis=0) 
            if np.isnan(q[0]): 
                q[0] = float('inf')
            D[m] = q 

        M = np.array(list(self.prefix2score_map.values()))
        idiff_mean,nvar_mean = np.mean(M,axis=0)
        return D,idiff_mean,nvar_mean 

    def fetch_scores_by_type(self,measure_type): 
        assert measure_type in DEFAULT_CYMETRIC_SUFFIX_MEASURE_TYPES
        D = dict() 
        for k,v in self.prefix2mcs_map.items(): 
            D[k] = v[measure_type]
        return D 

    def process_prefix(self,prefix): 
        P = {} 
        for l in DEFAULT_CYMETRIC_SUFFIX_MEASURE_TYPES: 
            m_score,ck_score = self.mcs_complexity_at_prefix(prefix,l)
            P[l] = (m_score,ck_score) 

        d0 = self.index_of_occurrence_difference_setsize_at_prefix(prefix)
        v0 = self.normalized_variance_of_suffix_frequency_at_prefix(prefix) 
        return P,d0,v0

    def index_of_occurrence_difference_setsize_at_prefix(self,prefix): 
        indices = sorted(self.M_.M[prefix].keys())
        V = stdop_vec(np.array(indices),sub,np.float32)
        l = len(set(V))
        if l == 0: return float('inf')
        return l 

    # NOTE: duplicate code? 
    def normalized_variance_of_suffix_frequency_at_prefix(self,prefix): 

        prefix_info = self.M_.key_to_seq_by_type(prefix,"valfreq") 
        frequencies = np.array([prefix_info[i] for i in range(1,len(prefix_info),2)])

        m = np.mean(frequencies)

        f2 = np.abs(frequencies - m) 
        mdiff = np.max(frequencies) - np.min(frequencies) 
        if mdiff == 0: return 0 

        return np.sum(f2) / (mdiff * len(f2))

    """
    [0] MCS k-complexity by averaging the four: 
    - bool + contiguous|best + median|most frequent
    [1] best cyclical kernel score 
    """
    def mcs_complexity_at_prefix(self,prefix,measure_type): 
        assert measure_type in DEFAULT_CYMETRIC_SUFFIX_MEASURE_TYPES

        prefix_info = self.M_.key_to_seq_by_type(prefix,measure_type) 
        ##print("PREFIX: ",prefix_info)
        ## print("MCS LEN: ",len(prefix_info))
        ms = MCSSearch(prefix_info,self.cast_func)
        ms.search() 

        m_score = self.mcs_complexity_at_prefix_(ms,prefix,measure_type)
        
        ck_score = ms.mean_cyclical_kernel_score()  
        return m_score,ck_score

    def mcs_complexity_at_prefix_(self,ms,prefix,measure_type): 
        assert type(ms) == MCSSearch 

        d0 = ms.default_kcomplexity(diff_type="bool",\
            diff_type2="best",basis="median")

        d1 = ms.default_kcomplexity(diff_type="bool",\
            diff_type2="contiguous",basis="median")

        d2 = ms.default_kcomplexity(diff_type="bool",\
            diff_type2="best",basis="most frequent")

        d3 = ms.default_kcomplexity(diff_type="bool",\
            diff_type2="contiguous",basis="most frequent")

        return (d0 + d1 + d2 + d3) / (4 * len(self.S)) 
