from intigers.mdr_v2 import * 
from intigers.extraneous import * 
from intigers.tvec import *
from mini_dm.ag_ext import *  
from mini_dm.ngram import *

class MultiMetric:

    def __init__(self,l):
        self.l = np.array(l)
        assert is_vector(self.l)
        self.ref_index = 0 

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

        agv2 = APRNGGaugeV2(f,(fmin,fmax),0.5) 

        for c in chunks:  
            cs = c.flatten()
            agv2.assign_cycle(cs)
            agv2.measure_cycle(len(cs),\
                term_func=lambda l,l2: type(l) != type(None),\
                auto_frange=not set_frange,auto_prange=True,\
                do_cycle_update=False) 

        # return the measurements from 
        return np.array(agv2.measurements)

    # TODO 
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
    def diff_measures(self):
        # get diff0,diff1 vectors 
            # pairwise diff
        v0 = stdop_vec(self.l,sub,np.float32)
            # change in pairwise diff
        v1 = stdop_vec(v0,sub,np.float32) 
        trel0 = to_trinary_relation_v2(v0,None,True,False)
        trel1 = to_trinary_relation_v2(v1,None,True,False)

        q = contiguous_repr__sequence(trel0) 
        q2 = contiguous_repr_sequence(trel1) 
        return len(q),len(q2) 

    def summarize(self):
        return -1 