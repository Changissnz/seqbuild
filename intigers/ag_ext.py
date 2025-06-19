from morebs2.aprng_gauge import *
from morebs2.matrix_methods import is_bounds_vector,is_proper_bounds_vector,is_2dmatrix  
from types import MethodType,FunctionType
from .seq_struct import * 
from .extraneous import * 

def absdiff_match_func(i,i2):
    mm = dict() 
    x = abs(i2 - i)
    q = min(x)
    indices = np.where(x == q)[0]
    return list(i2[indices])

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

    @staticmethod 
    def measure_match(match_map):
        # measure the tie-size
        q = sum([len(v) - 1 for v in match_map.values()])

        # measure the most frequent intersection
        mmf = MinMaxFreq(match_map,False)
        while not mmf.fin: 
            mmf.count_one() 
        mmf.finalize_count()
        
        fx = mmf.nth_most_frequent(0)
        return q,fx 

    @staticmethod
    def match_two_intseq(i1,i2,match_func): 
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
    The `seg_length` is assigned the value of the mean pairwise 
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
    
"""
Measures categorical entropy over rows or columns of matrix `m`. The variable 
`sl_info` is for the `seg_length` parameter of <APRNGGaugeV2.std_cat_entropy>. 
It is one of five types: None, int (denumerator for max pairwise distance), 
vector<int> (denumerator for each row or column), float (literal), and vector<float>. 
"""
def APRNGGaugeV2__matrix_cat_entropy(m,franges,is_rowwise:bool=True,is_local_frange:bool=True,\
    sl_info = None,count_type="absdiff",round_depth:int=5):
    assert is_2dmatrix(m)
    m_ = m.T if not is_rowwise else m 
    f = None 
    i = None  

    if type(franges) == type(None): 
        if is_local_frange: 
            franges = []
            for x in m_: 
                franges.append((np.min(x),np.max(x))) 
            franges = np.array(franges,dtype=np.int32) 
        else: 
            franges = (np.min(m_),np.max(m_))

    if is_bounds_vector(franges): 
        assert is_proper_bounds_vector(franges)
        assert franges.shape[0] == len(m_)
        f,i = franges[0],0 
    else: 
        assert type(franges) == tuple
        assert len(franges) == 2 
        f = franges 

    j = 0 if is_vector(sl_info) else None
    sl = sl_info  

    ag = APRNGGaugeV2(None,f,pradius=5)
    lx = [] 
    for x in m_: 
        iseq = IntSeq(x) 

        if type(i) != type(None): 
            f = franges[i] 
            ag.reload_var("frange",tuple(f) )

        if type(j) != type(None):
            sl = sl_info[j] 
            assert sl > 0

        sl_ = sl 
        if type(sl_) in {int,np.int32,np.int64}:
            sl_ = (f[1] - f[0]) / sl 

        ce = ag.std_cat_entropy(iseq,seg_length=sl_,start_value=f[0],\
            count_type=count_type)
        lx.append(ce) 

        if type(i) != type(None): 
            i += 1 

        if type(j) != type(None): j += 1 

    return np.round(np.array(lx),round_depth)
