from morebs2.aprng_gauge import *
from morebs2.matrix_methods import is_bounds_vector,is_proper_bounds_vector,is_2dmatrix  
from types import MethodType,FunctionType
from .minmax_freq import *  
from .iseq import * 
from intigers.process_seq import stdop_vec
from intigers.extraneous import zero_div0  

"""
chooses the element of i2 closest in absolute distance 
to float i. 

i := float 
i2 := vector 
"""
def absdiff_match_func(i,i2):
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

    def update_cycle(self,max_size:int,term_func):
        self.cycle = self.cycle_one(max_size,term_func)
        return

    def measure_cycle(self,max_size,\
        term_func=lambda l,l2: type(l) == type(None),\
        auto_frange:bool=False,auto_prange:bool=False,\
        do_cycle_update:bool=False):

        if do_cycle_update:
            self.update_cycle(max_size,term_func) 

        if auto_prange:
            sz0,sz1 = min(self.cycle),max(self.cycle)
            self.pradius = APRNGGaugeV2.default_pradius(sz0,sz1,len(self.cycle)) 

        q = super().measure_cycle(max_size,term_func,auto_frange) 
        return q

    # CAUTION: does not check if dim. (d0,d1) is possible. 
    def output_to_matrix(self,d0,d1):
        m = []
        term_func=lambda l,l2: type(l) != type(None)

        for _ in range(d0):
            self.update_cycle(d1,term_func)
            assert len(self.cycle) == d1  
            m.append(deepcopy(self.cycle))  
        return np.array(m) 

    """
    measures a matrix `m` for qualities. Stores 
    values into dictionary 

    D: axis -> index -> 
        [0] (coverage,unidirectional weighted point distance)
        [1] entropy::float 
        [2] (min pairwise difference, max pairwise difference,mean pairwise difference)
    """
    def measure_matrix_chunk(self,m,d0,d1,axes={0,1}): 
        assert axes.issubset({0,1}) 

        if not is_2dmatrix(m): 
            m = self.output_to_matrix(d0,d1)
        else: 
            assert m.shape == (d0,d1)  

        # do along each axis, measure cycle 
        D = dict() 
        while len(axes) > 0: 
            A = axes.pop() 
            D = APRNGGaugeV2.measure_matrix__rc_wise(self,D,m,A)
        return D 

    """
    auxiliary method used by method<measure_matrix_chunk> to 
    calculate values for qualities of matrix `m`. 

    D := dict, stores values, output from this method. 
    m := matrix of values
    a := axis, either 0 or 1. 
    """
    @staticmethod 
    def measure_matrix__rc_wise(AG,D,m,a): 
        assert a in {0,1}

        D[a] = dict()
        d = m.shape[a]

        def next_vec(i):
            q = m[:,i] if a == 1 else m[i,:] 
            return q 

        tfunc = lambda l,l2: type(l) != type(None)
        for i in range(d): 
            q = next_vec(i) 
            iseq = IntSeq(q)
            ent0 = AG.std_cat_entropy(iseq,seg_length=None,\
            start_value=None,count_type="absdiff")
            pdiff0 = APRNGGaugeV2.pairwise_diff_metrics(iseq) 

            AG.assign_cycle(q) 
            m0 = AG.measure_cycle(len(q),\
                term_func=tfunc,auto_frange=True,\
                auto_prange=True,do_cycle_update=False)
            D[a][i] = (m0,ent0,pdiff0)
        return D 

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
        #assert type(i1) == IntSeq
        return i1.match_map(i2,match_func) 

    """
    return: 
    - min,max,mean of pairwise differences from is1
    """
    @staticmethod
    def pairwise_diff_metrics(is1):
        #assert type(is1) == IntSeq
        if len(is1) < 2:
            return None,None,None 
        q = uwpd(np.array(is1.l),pairwise_op=lambda x1,x2: np.abs(x2 - x1),accum_op=None)
        return min(q),max(q), sum(q) / len(q) 

    @staticmethod 
    def default_pradius(min_value,max_value,sample_size):
        assert sample_size >= 1 
        assert max_value >= min_value
        return (max_value - min_value) / sample_size

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
            _,_,seg_length = APRNGGaugeV2.pairwise_diff_metrics(is1)

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