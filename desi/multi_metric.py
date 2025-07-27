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
    def vector_measures(self):
        # get diff0,diff1 vectors 
        v0 = stdop_vec(self.l,sub,np.float32)
        v1 = stdop_vec(v0,sub,np.float32) 
        return (v0,v1) 

    def summarize(self):
        return -1 