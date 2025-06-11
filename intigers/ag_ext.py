from morebs2.aprng_gauge import *
from .seq_struct import * 

def absdiff_match_func(i,i2):
    mm = dict() 
    x = abs(i2 - i)
    q = min(x)
    indices = np.where(x.l == q)[0]
    return i2[indices]

# TODO: test 
class APRNGGaugeV2(APRNGGauge):

    def __init__(self,aprng,frange,pradius:float):
        super().__init__(aprng,frang,pradius)

    def measure_cycle(self,max_size,\
        term_func=lambda l,l2: type(l) == type(None)): 
        q = super().measure_cycle(max_size,term_func)
        return q

    def measure_match(self,match_map):
        # measure the tie-size
        q = sum([len(v) - 1 for v in match_map.values()])

        # measure the most frequent intersection
        mmf = MinMaxFreq(match_map,False)
        mmf.count_one() 
        mmf.finalize_count()
        
        q = mmf.sorted_counts[-1][1]
        fx = []

        for x in mmf.sorted_counts[::-1]:
            if x[1] != q: break 
            fx.append(x[0]) 
        return q,fx 

    def match_two_intseq(self,i1,i2,match_func): 
        assert type(i1) == type(i2)
        assert type(i1) == IntSeq
        return i1.match_map(i2,match_func) 

    """
    return: 
    - min,max,mean of pairwise differences from is1
    """
    @staticmethod
    def pairwise_diffs(is1):
        assert type(is1) == IntSeq

        if len(is1) < 2:
            return None,None,None 

        c = 0
        dx = 0 
        minnie,maxie = float('inf'),-float('inf')
        for i in range(len(is1)-1):
            q = is1[i]
            for j in range(i+1,len(is1)):
                x = abs(q - is1[j])
                c += x
                if x < minnie: minnie = x
                if x > maxie: maxie = x
                dx += 1
        return minnie,maxie, c / dx

    """
    standard categorical entropy measures the number of 
    contiguous pairwise difference in category. Categories 
    are assigned by using the `start_value` (the minumum 
    value of `is1` if `None`) as the 0 category. An element 
    x is labeled the category `(x - start_value) / seg_length`. 
    The `seg_length` is assigned the value of the minimum pairwise 
    difference of `is1`. 
    """
    def std_cat_entropy(self,is1,seg_length:float=None,start_value=None,\
        count_type="absdiff"):
        assert count_type in {"absdiff","equals"}

        if len(is1) < 2: 
            return 0.0 

        if type(seg_length) == type(None):
            seg_length,_,_ = APRNGGaugeV2.pairwise_diffs(is1)

        if seg_length == 0.0: return 0.0 

        catvec = is1.diffcat_vec(seg_length,start_value)
        if len(np.unique(catvec)) < 2: return 0.0 

        c = 0
        for i in range(len(catvec) - 1):
            q = abs(catvec[i] - catvec[i+1])
            if q == 0.0: continue 

            c += 1 if count_type == "equals" else q 
        
        qc = max(catvec) if count_type == "absdiff" else 1 
        return c / (qc * (len(catvec) - 1)) 