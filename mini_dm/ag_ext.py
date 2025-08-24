from morebs2.aprng_gauge import *
from morebs2.matrix_methods import is_bounds_vector,is_proper_bounds_vector,is_2dmatrix  
from types import MethodType,FunctionType
from .minmax_freq import *  
from .iseq import * 
from intigers.process_seq import stdop_vec
from intigers.extraneous import zero_div0,to_trinary_relation_v2
from intigers.intfactor import * 
from intigers.poly_output_fitter_ import * 

#------------------------ auxiliary functions used by <APRNGGaugeV2> for 
#------------------------ measures such as matching and summarization 

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

def summarize_matrix(m):
    q = np.array([np.min(m,axis=0),\
        np.max(m,axis=0),np.mean(m,axis=0),\
        np.var(m,axis=0)]) 
    return q 

#---------------- methods used for sequence modification according 
#---------------- to metrics such as coverage and u.w.p.d. 

def trinary_vec_for_element(S,i):
    dv = to_trinary_relation_v2(np.array(S),S[i],True,False) 
    return dv 

"""
calculates a sequence D, named ranged-delta decomposition, 
w.r.t. the i'th element of `S`. The direction `d` of changing 
S[i] is either -1 or 1. The elements of D  are of the form: 
    ((start additive,end additive), rate of change). 
The super-bound of these additives w.r.t. S[i] is 
    [rv[0] - S[i],0] if d == -1 and 
    [0,rv[1] - S[i]] if d == 1. 

The sequence D is a partitioning used to measure the 
change of unidirectional weighted point distance 
from the change of S[i] over spans of additives to 
S[i]. 
"""
def ranged_delta_decomposition(S,i,rv=None,d=1):
    assert d in {-1,1} 

    mn,mx = np.min(S),np.max(S) 
    if type(rv) == type(None):
        rv = (mn,mx) 
    assert is_valid_range(rv,False,True) 
    assert rv[0] <= mn <= rv[1]
    assert rv[0] <= mx <= rv[1]

    dvec = trinary_vec_for_element(S,i) 

    iset0,iset1 = np.where(dvec == -1)[0],\
        np.where(dvec == 1)[0]

    # d = -1 -> (-1,-1),(1,1) 
    # d = 1 -> (-1,1),(1,-1) 
    isets = [iset0,iset1]

    dx = d if d == 1 else 0 
    max_diff = rv[dx] - S[i]

    indices = isets[dx]
    diff_vec = [] 
    for j in indices:
        diff = abs(S[j] - S[i])
        diff_vec.append(diff)
    diff_vec = sorted(diff_vec) 

    start_index = 0.0 
    iset_count = [len(isets[0]),len(isets[1])] 
    iset_ = set(np.where(dvec == 0)[0]) - {i}
    iset_count[(dx + 1) % 2] = iset_count[(dx + 1) % 2] + \
        len(iset_)

    decomp = []
    for dv in diff_vec:
        end_index = round(dv * d,5)
        r = (start_index,end_index) 
        per = iset_count[(dx + 1) % 2] - iset_count[dx]
        decomp.append((r,per)) 
        start_index = end_index 

        if iset_count[dx] == 0: 
            continue
        iset_count[dx] -= 1 
        iset_count[(dx+1)%2] += 1 
   
    if start_index != max_diff:
        r = (start_index,max_diff) 
        per = iset_count[(dx + 1) % 2] - iset_count[dx]
        decomp.append((r,per))

    return decomp 

"""
determines the index of the delta range in `rdd`, 
a ranged-delta decomposition, for the wanted change 
`delta` to occur. If there is no possibility for `delta` 
to occur, given `rdd`, then the algorithm outputs -1 
as the index. 

return: 
- index of `rdd` 
- range of change for `rdd`[index]
- remaining of `delta` at that index in `rdd` 
"""
def rdd_index_for_delta(rdd,delta): 
    i = -1
    s = 0
    dstat = delta > 0 
    for (j,x) in enumerate(rdd):
        d = abs(x[0][1] - x[0][0]) 

        if x[1] == 0: continue 
        dx = d * x[1] 

        s1 = s + dx
        s_ = sorted([s,s1])

        sdiff = s1 - s 
        if dstat: 
            if delta <= sdiff:
                return j,[s,s1],delta
        else:
            if delta >= sdiff: 
                return j,[s,s1],delta 

        delta = delta - dx 
        s = s1  
    return i,None,None 


"""
adjusts the i'th value of S so that the unidirectional weighted point 
distance (u.w.p.d) after the adjustment is of difference `c` to the original 
u.w.p.d. Depending on the magnitude and sign of `c`, it may not be possible 
to achieve that change in u.w.p.d. In these cases, method outputs None. 
The new value of S[i], if `c` is possible to achieve, is in the confines of 
the super-range `rv`, which is (min(S),max(S)) if set to None. 
"""
def adjust_for_uwpd_change(S,i,c,rv=None,d_priority=1,recurse:bool=True): 
    rdd = ranged_delta_decomposition(S,i,rv,d_priority)

    di,rx,delta = rdd_index_for_delta(rdd,c) 

    if di == -1:
        if not recurse: 
            return None
        return adjust_for_uwpd_change(S,i,c,rv,d_priority * -1,False) 

    rx_ = rdd[di][0]
    ratio = delta / (rx[1] - rx[0])

    S_ = np.asarray(S,dtype=float) 

    q = delta / rdd[di][1]
    if rx_[1] < rx_[0]:
        q = -q 

    S_[i] = S_[i] + rx_[0] + q
    return S_ 

#----------------------------------- extension of <morebs2.APRNGGauge> 

CHAINTEST_DESCRIPTOR = "mean(min),mean(max)\nvar(min),var(max),\nmean(F_i - F_(i+1)),mean(F_(i+1) - F_i),\nvar(F_i - F_(i+1)),var(F_(i+1) - F_i)"

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
        self.clear_cycle() 
        q = self.cycle_one(max_size,term_func)
        self.cycle = q 
        return

    """
    NOTE: 
    in some cases with excessively large values, 
    such as with np.int32 values that are almost the 
    maximum np.int32, the output measurement 
        [0] coverage
        [1] normalized unidirectional weighted point distance
    may not make sense, especially [1] when [1] is not in the 
    range [0.,1.] (that is the appropriate range for normalized 
    u.w.p.d.).

    When the u.w.p.d. is erroneous, implying the involvement of 
    an excessively large value, there is usually a warning during 
    Python program execution that reads 
        >>> RuntimeWarning: overflow encountered
    """
    def measure_cycle(self,max_size,\
        term_func=lambda l,l2: type(l) == type(None),\
        auto_frange:bool=False,auto_prange:bool=False,\
        do_cycle_update:bool=False):

        if do_cycle_update:
            self.update_cycle(max_size,term_func) 

        if auto_prange:
            sz0,sz1 = min(self.cycle),max(self.cycle)
            self.pradius = APRNGGaugeV2.default_pradius(sz0,sz1,len(self.cycle))
            self.pradius = max([self.pradius,10 ** -5])

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
    measures a matrix `m` along its rows or columns 
    for qualities. Stores values into dictionary 

    D: axis -> index -> 
        [0] (coverage,unidirectional weighted point distance)
        [1] entropy::float 
        [2] (min,max,mean) pairwise difference
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
            D = APRNGGaugeV2.measure_matrix__rc_wise(self,D,m,A,\
                auto_frange=True,auto_prange=True)
        return D 

    """
    auxiliary method used by method<measure_matrix_chunk> to 
    calculate values for qualities of matrix `m`. 

    D := dict, stores values, output from this method. 
    m := matrix of values
    a := axis, either 0 or 1. 
    """
    @staticmethod 
    def measure_matrix__rc_wise(AG,D,m,a,\
        auto_frange=True,auto_prange=True): 
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
                term_func=tfunc,auto_frange=auto_frange,\
                auto_prange=auto_prange,do_cycle_update=False)
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
        assert issubclass(type(i1),IntSeq)
        assert issubclass(type(i2),IntSeq)
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
        seg_length = abs(seg_length)

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
    lightweight alternative measurement to the others. 

    The chain test measures the continuity of qualities 
    shared from one contiguous chunk to the other. There 
    are a total of `num_iter` elements, each chunk having 
    `chain_length` elements (or remainder). 

    For chunks C_1,C_2,...C_j, calculates the mean and 
    variance of these qualities: 
    - minumum
    - maximum 

    The factor-frequency difference maps (see <ISFactorSetOps>) 
    from the chunks C_1,C2,...,C_j are 
    D_1,...,D_(j-1). 

    And their compressed float forms are 
    f_1,...,f_j (see method<accum_factor_frequency_map>).

    Output from the test is then a 4 x 2 matrix with values 
    for: 

    mean(min),mean(max),
    var(min),var(max),
    mean(D_i - D_(i+1)),mean(D_(i+1)-D_i),   
    var(D_i - D_(i+1)),var(D_(i+1)-D_i).
    """
    def chaintest(self,num_iter,chain_length):
        assert type(num_iter) in {int,np.int32,np.int64} 
        assert type(chain_length) in {int,np.int32,np.int64} 

        assert 2 <= chain_length < num_iter 

        num_chunks = ceil(num_iter / chain_length)
        tx = [] # get var,mean at end 
        tx2 = [] 
        
        def factor_func(cl):
            X = np.array([self.aprng() for _ in range(cl)]) 
            X = [modulo_in_range(x,[-NPINT32_MAX+1,NPINT32_MAX]) for x in X]
            X = np.array(X,dtype=np.int32)


            m0,m1 = np.min(X),np.max(X) 
            tx.append((m0,m1)) 
            isfso = ISFactorSetOps(X)
            isfso.factor_count_()
            return isfso 

        isf = factor_func(chain_length)
        num_chunks -= 1 
        for i in range(num_chunks): 
            l = min([num_iter,chain_length]) 
            num_iter -= l 

            isf2 = factor_func(l)
            dx = isf - isf2 
            dx2 = isf2 - isf 

            n0 = accum_factor_frequency_map(dx)
            n1 = accum_factor_frequency_map(dx2) 
            tx2.append((n0,n1)) 

        tx,tx2 = np.array(tx),np.array(tx2)

        minmax_mean = np.mean(tx,axis=0)
        minmax_var = np.var(tx,axis=0)
        fdiff_mean = np.mean(tx2,axis=0) 
        fdiff_var = np.var(tx2,axis=0) 

        return np.array([minmax_mean,\
            minmax_var,fdiff_mean,fdiff_var])