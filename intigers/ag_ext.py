from morebs2.aprng_gauge import *
from types import MethodType,FunctionType
from .seq_struct import * 
from .extraneous import * 

def absdiff_match_func(i,i2):
    mm = dict() 
    x = abs(i2 - i)
    q = min(x)
    indices = np.where(x.l == q)[0]
    return i2[indices]

# TODO: test 
class APRNGGaugeV2(APRNGGauge):

    def __init__(self,aprng,frange,pradius:float):
        super().__init__(aprng,frange,pradius)
        self.catvec = None 

    def reload_var(self,varname,varvalue):
        assert varname in {"aprng","frange","pradius"}

        if varname == "aprng":
            assert type(varvalue) in {MethodType,FunctionType}
        elif varname == "frange":
            assert type(varvalue) in {tuple,list}
            assert len(varvalue) == 2 
            assert varvalue[0] <= varvalue[1]
        else: 
            assert type(varvalue) == float 
            assert varvalue >= 0.0 and varvalue <= 1.0 
        
        setattr(self,varname,varvalue)

    def clear_cycle(self):
        self.cycle = None 

    def measure_cycle(self,max_size,\
        term_func=lambda l,l2: type(l) == type(None),\
        auto_frange:bool=False): 
        q = super().measure_cycle(max_size,term_func,auto_frange) 
        return q

    def measure_match(self,match_map):
        # measure the tie-size
        q = sum([len(v) - 1 for v in match_map.values()])

        # measure the most frequent intersection
        mmf = MinMaxFreq(match_map,False)
        mmf.count_one() 
        mmf.finalize_count()
        fx = mmf.nth_most_frequent(0)
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
        q = uwpd(np.array(is1.l),pairwise_op=lambda x1,x2: np.abs(x2 - x1),accum_op=None)
        return min(q),max(q), sum(q) / len(q) 

    """
    standard categorical entropy measures the number of 
    contiguous pairwise differences in category. Categories 
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
            _,_,seg_length = APRNGGaugeV2.pairwise_diffs(is1)

        if seg_length == 0.0: return 0.0 
        self.catvec = is1.diffcat_vec(seg_length,start_value)
        if len(np.unique(self.catvec)) < 2: return 0.0 

        c = 0
        for i in range(len(self.catvec) - 1):
            q = abs(self.catvec[i] - self.catvec[i+1])
            if q == 0.0: continue 

            c += 1 if count_type == "equals" else q 
        
        qc = max(self.catvec) if count_type == "absdiff" else 1 
        return c / (qc * (len(self.catvec) - 1)) 

    def cycle_multvec(self):
        return stdop_vec(self.cycle,zero_div0,cast_type=np.float32)