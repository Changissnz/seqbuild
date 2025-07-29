from intigers.mdr_v2 import * 
from intigers.extraneous import * 
from intigers.tvec import *
from intigers.intfactor import * 
from mini_dm.ag_ext import *  
from mini_dm.ngram import *
from morebs2.seq_repr import * 

"""
calculates the Kolmogorov complexity via 
size of the contiguous representation of the 
1st-order difference and 2nd-order difference; 
these differences are in trinary form before 
being mapped into contiguous representation.

Outputs a pair of non-negative integers 
(p1,p2), in which p1 is the size of the contiguous 
representation of the 1st-order difference (in trinary form) 
and p2 likewise for the 2nd-order difference. 
"""
def trinary_diffvec_kcomplexity(l):
    # get diff0,diff1 vectors 
        # pairwise diff
    v0 = stdop_vec(l,sub,np.float32)
        # change in pairwise diff
    v1 = stdop_vec(v0,sub,np.float32) 
    trel0 = to_trinary_relation_v2(v0,None,True,False)
    trel1 = to_trinary_relation_v2(v1,None,True,False)

    q = contiguous_repr__sequence(trel0) 
    q2 = contiguous_repr__sequence(trel1) 

    return len(q),len(q2) 

"""
Modular characteristic of integer sequence `l`. 
Calculates a map consisting of keys that are 
the |`l`|'th most frequent factors of the elements 
of `l`, and corresponding value the Kolmogorov 
complexity of contiguous condensed representations 
for 1st and 2nd order differences. 
"""
class ModChar: 

    def __init__(self,l): 
        assert type(l) == list or is_vector(l) 
        self.l = np.array(l,dtype="int32") 
        self.factors,self.complexity_map = None,None 
        self.init_factors()
        return

    def init_factors(self): 
        isfso = ISFactorSetOps(self.l,int_limit=DEFAULT_INT_MAX_THRESHOLD)
        isfso.factor_count_() 
        dx2 = isfso.dsort(pkeys=None)[::-1] 

        l = len(self.l)
        dx_ = [dx2_[0] for dx2_ in dx2[:l]] 
        self.factors = dx_ 
        if 2 not in self.factors: self.factors.insert(0,2) 

    """
    main method 
    """
    def modular_complexity(self): 
        self.complexity_map = dict()
        for f in self.factors: 
            l2 = np.array([l_ % f for l_ in self.l])
            self.complexity_map[f] = trinary_diffvec_kcomplexity(l2)

"""
similar to a swiss-knife but for measuring qualities 
of sequence `l`. 
"""
class MultiMetric:

    def __init__(self,l):
        self.l = np.array(l,dtype=np.float64)
        assert is_vector(self.l)
        self.modcomplex_map = None 

    def load_mc_map(self):
        mc = ModCond(self.l) 
        mc.modular_complexity()
        self.modcomplex_map = mc.complexity_map 

    """
    produces sequence of two-dimensional measures 
    (coverage, normalized unidirectional weighted point distance) 
    for subsequence S of `l`. S starts at `ref_index` of `l` and 
    spans for `length` number of elements in increasing index order. The 
    positive integer `l2` is the subsequence length for the n-gram 
    iteration over S.
    """
    def agv2_measures__ngrammer(self,ref_index,length,\
            l2,set_frange:bool=True): 
        assert l2 <= length 

        # init the <NGrammer> 
        sv = subvec(self.l,ref_index,length) 

        ng = NGrammer2D(sv,l2,(1,l2))
        chunks = ng.one_cycle(ref_index=0)

        def f(): return None

        fmin,fmax = 0.,1.
        if set_frange: 
            fmin,fmax = np.min(self.l),np.max(self.l) 
        
        start_value = np.min(self.l)

        agv2 = APRNGGaugeV2(f,(fmin,fmax),0.5) 

        ce_vec = []
        for c in chunks:  
            cs = c.flatten()
            agv2.assign_cycle(cs)
            agv2.measure_cycle(len(cs),\
                term_func=lambda l,l2: type(l) != type(None),\
                auto_frange=not set_frange,auto_prange=True,\
                do_cycle_update=False) 

            is1 = IntSeq(cs) 
            ce = agv2.std_cat_entropy(is1,seg_length=None,start_value=start_value,\
                count_type="absdiff")
            ce_vec.append(ce) 

        # return the measurements from 
        q = np.array(agv2.measurements)
        ce_vec = np.array(ce_vec)
        ce_vec = ce_vec.reshape((len(ce_vec),1))
        return np.hstack((q,ce_vec))

    """
    see method<trinary_diffvec_kcomplexity>
    """
    def diff_measures(self):
        return trinary_diffvec_kcomplexity(self.l) 

    """
    Kolmogorov complexity by calculation involving most common 
    subsequence (actually, the technical term is subsequence of 
    median frequency).  

    Value falls in range [0.,1.] and is directly proportional to 
    the complexity (measured by the error of the subsequence fitted 
    to the target sequence `l`) of representing `l`. 
    """
    def mcs_kcomplexity(self):
        d = MCS_kcomplexity(self.l,float,diff_type="bool",\
            diff_type2="best",basis="median")
        return zero_div0(d,len(self.l))

    """
    Provides a summarization of sequence `l`. The `ngram` argument 
    specifies the length of contiguous subsequences from `l` that are 
    measured for particular qualities (see description of these values 
    below `return`). 

    return: 
    [0] (coverage, normalized unidirectional weighted point distance, 
        standard categorical entropy)
    [1] Kolmogorov complexity of representing 1st-order difference,
        Kolmogorov complexity of representing 2nd-order difference,
    [2] Kolmogorov complexity of difference b/t `l` and its most 
        common subsequence. 
    """
    def summarize(self,ngram):
        m1 = self.agv2_measures__ngrammer(0,len(self.l),\
            ngram,set_frange=True)
        m1 = np.mean(m1,axis=0)
        m2 = self.diff_measures()
        m3 = self.mcs_kcomplexity()
        return m1,m2,m3