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

NOTE: modular characteristic always includes 2, so the 
      number of factors may be either |`l`| or |`1`| + 1. 
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

    """
    loads map describing modular characteristic of `l`. 
    """
    def load_mc_map(self):
        mc = ModChar(self.l) 
        mc.modular_complexity()
        self.modcomplex_map = mc.complexity_map 

    """
    produces sequence of three-dimensional measures 
    (coverage, normalized unidirectional weighted point distance,
    categorical entropy), 
        (float,float,float)
    for subsequence S of `l`. S starts at `ref_index` of `l` and 
    spans for `length` number of elements in increasing index order. The 
    positive integer `l2` is the subsequence length for the n-gram 
    iteration over S.
    """
    def agv2_measures__ngrammer(self,ref_index,length,\
            l2,set_frange:bool=True,external_frange=None): 
        assert l2 <= length, "length,l2 is {}".format((length,l2))

        # init the <NGrammer> 
        sv = subvec(self.l,ref_index,length) 

        ng = NGrammer2D(sv,l2,(1,l2))
        chunks = ng.one_cycle(ref_index=0)

        def f(): return None

        fmin,fmax = 0.,1.
        if set_frange: 
            fmin,fmax = np.min(self.l),np.max(self.l) 
        
        if type(external_frange) != type(None):         
            assert is_valid_range(external_frange,False,False)
            fmin,fmax = external_frange

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
    subsequence (actually, the technical wording is subsequence of 
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

        * every value in [0]-[2] is a float. 
    """
    def summarize(self,ngram,condense_ngram_output:bool=True,set_frange:bool=True):
        m1 = self.agv2_measures__ngrammer(0,len(self.l),\
            ngram,set_frange=set_frange)
        if condense_ngram_output:
            m1 = np.mean(m1,axis=0)
            
        m2 = self.diff_measures()
        m3 = self.mcs_kcomplexity()
        return m1,m2,m3


#------------------------------------------------------------------------------------

# below methods are used for interface command QUALTEST. 

"""
method used to average the sequence of output elements 
from method<gauge_generator__MultiMetric>.
"""
def average_MultiMetric_summaries(R):
    l = len(R) 
    if l == 0: return None 

    q = R.pop(0) 
    N0,N1,N2 = q[0],np.array(q[1]),q[2] 
    while len(R) > 0: 
        q = R.pop(0) 
        n0,n1,n2 = q[0],np.array(q[1]),q[2] 
        N0 += n0 
        N1 += n1 
        N2 += n2 
    
    return N0 / l, N1 / l, N2 / l 

"""
a heavyweight function used to measure output qualities of a 
pseudo-random number generator `prg`. 

NOTE: function conducts intensive calculations, by way of n-gram 
iterations. On personal computing systems, the argument `num_iter` 
normally should not exceed 1000. This rule translates to this method 
having a size limit for measuring sequences. 

Function outputs the following measures. 
[0] average or list of measures from method<MultiMetric.summarize>,
[1] modular characteristic map,
[2] element frequency map.
"""
def gauge_generator__MultiMetric(prg,num_iter,gauge_depth:int,deg_vec:list,\
    set_frange:bool=True,condense_ngram_output:bool=True): 
    assert num_iter >= 5 
    
    stat0 = type(gauge_depth) == type(None) 
    stat1 = type(deg_vec) == type(None) 

    assert stat0 or stat1 and not (stat0 and stat1) 

    if stat1: 
        base_ngram = 10 
        if num_iter < base_ngram: 
            base_ngram = int(round(num_iter/2)) 
    else: 
        base_ngram = None 

    index = 0 if stat0 else None 

    def next_ngram(ngram_,index_): 
        if stat1: 
            ngram_ += 1 
            if ngram_ >= num_iter:
                return None,None

            if ngram_ > base_ngram + gauge_depth: 
                return None,None 

            return ngram_,None

        if index_ >= len(deg_vec): 
            return None,None
        v = deg_vec[index_]
        index_ += 1 
        return v,index_

    # collect all elements into list 
    q = [] 
    for _ in range(num_iter): 
        q.append(prg())

    mm = MultiMetric(q)
    stat = True 
    
    ngx = base_ngram 
    R = [] 
    while stat: 
        ngx,index = next_ngram(ngx,index) 
        stat = not type(ngx) == type(None) 
        if not stat: continue 
        A = mm.summarize(ngx,condense_ngram_output=True,set_frange=set_frange)
        R.append(A) 

    fm = vec_to_frequency_map(np.array(q,dtype=int))

    cmap = None 
    try: 
        mm.load_mc_map() 
        cmap = mm.modcomplex_map
    except: 
        pass 

    if condense_ngram_output: 
        R = average_MultiMetric_summaries(R)

    return R,cmap,fm 

"""
a comparator between generators `prg` and `prg2`. Function is 
based on function<gauge_generator__MultiMetric>. Calculates 
the difference 

`gauge_generator__MultiMetric(prg,...) - gauge_generator__MultiMetric(prg2,...)`. 
"""
def cmp_generators__MultiMetric(prg,prg2,num_iter,gauge_depth,deg_vec,\
    set_frange:bool=True): 
    q0 = gauge_generator__MultiMetric(prg,num_iter,gauge_depth,deg_vec,set_frange)
    q1 = gauge_generator__MultiMetric(prg2,num_iter,gauge_depth,deg_vec,set_frange)

    def diff_mm_output():
        x0,x1 = q0[0],q1[0] 

        d0 = list(x0[0]) 
        d0.extend(x0[1]) 
        d0.append(x0[2]) 
        d0 = np.array(d0) 

        d1 = list(x1[0]) 
        d1.extend(x1[1]) 
        d1.append(x1[2]) 
        d1 = np.array(d1) 
        return d0 - d1 

    def diff_mc_map():
        x0,x1 = q0[1],q1[1] 
        if type(x0) == type(None) or type(x1) == type(None): 
            return dict() 

        ks = set(x0.keys()) | set(x1.keys()) 

        diff_map = dict() 

        for k in ks: 
            stat0 = k in x0 
            stat1 = k in x1 

            if stat0 and stat1: 
                d0,d1 = x0[k],x1[k] 
                diff_map[k] = tuple(np.abs(np.array(d0)-np.array(d1)))
            elif stat0:
                diff_map[k] = x0[k] 
            elif stat1:
                diff_map[k] = x1[k] 

        return diff_map 

    def diff_fmap(): 
        x0,x1 = q0[2],q1[2] 
        
        l0,l1 = len(x0),len(x1) 
        dx0 = l0 - l1 

        fx = set(x0.keys()) | set(x1.keys()) 

        dx1 = dict()
        for fx_ in fx: 
            c0 = x0[fx_] if fx_ in x0 else 0 
            c1 = x1[fx_] if fx_ in x1 else 0 
            dx1[fx_] = abs(c0-c1) 

        return dx1 

    dx0 = diff_mm_output()
    dx1 = diff_mc_map() 
    dx2 = diff_fmap() 
    return dx0,dx1,dx2 