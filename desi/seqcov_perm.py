from mini_dm.ag_ext import * 
from morebs2.numerical_generator import prg_partition_for_float
from .multi_metric import * 

"""
super-class for <SeqUWPDPermuter> and <SeqCoveragePermuter>. 
Class is used to apply changes to a numerical sequence for 
particular changes in the sequence's qualities (specifically 
coverage, unidirectional weighted point distance (uwpd), and 
modular uwpd). The `c_delta` value is the magnitude of change to 
be applied to `sequence`. 
"""
class AGV2SeqQualPermuter:

    def __init__(self,sequence,c_delta,super_range,prg): 
        assert is_vector(sequence) 
        assert type(c_delta) in {float,np.float32,np.float64}
        assert -1. <= c_delta <= 1. 
        assert is_valid_range(super_range,False,False) or type(super_range) == type(None)
        assert type(prg) in {MethodType,FunctionType} 

        self.l = sequence 
        self.c_delta = c_delta
        self.srange = super_range
        self.prg = prg 
        return 

"""
permuter that changes a base `sequence` according to the wanted 
change in unidirectional weighted point distance, `uwpd_delta`. 
The `uwpd_delta` value is a ratio (in [0.,1.]). 

The `super_range` is an ordered 2-tuple, and all values from `sequence` 
should be contained in it. This 2-tuple is used to calculate maximum 
possible uwpd. 

NOTE: 
If `modulus` is not set to `None`, instead a float f_i, then class instance 
is used to permute the `sequence` by its modular form, `sequence` % f_i. 
This implies the class serves dual purposes of permuting sequences by 
their absolute or modular unidirectional weighted point distance. 
"""
class SeqUWPDPermuter(AGV2SeqQualPermuter):

    def __init__(self,sequence,uwpd_delta,super_range,prg,modulus=None): 
        super().__init__(sequence,uwpd_delta,super_range,prg)
        self.modulus = modulus 
        self.l2 = None 
        self.l_ = None 
        self.preproc() 
        self.null_ctr = 0 
        return

    def preproc(self):
        # get max uwpd
        mxx = max_float_uwpd(len(self.l),self.srange)
        self.change_balance = mxx * self.c_delta
        
        self.l_ = deepcopy(self.l)         
        if type(self.modulus) != type(None):
            self.l = [l_ % self.modulus for l_ in self.l]   
        self.l2 = deepcopy(self.l)

    def modular_vec_to_closest_vec(self,mod_vec):
        cvec = []
        for (i,l_) in enumerate(self.l_): 
            q = l_ // self.modulus
            cvec.append(q * self.modulus + mod_vec[i])
        return np.array(cvec) 

    """
    main method 
    """
    def apply(self,null_limit=10): 
        if self.change_balance > 0.:
            def fx():
                return self.change_balance > 0.
        else: 
            def fx(): 
                return self.change_balance < 0. 

        while fx() and self.null_ctr < null_limit: 
            self.__next__()

        if type(self.modulus) == type(None): 
            return self.l2 
        return self.modular_vec_to_closest_vec(self.l2)

    def __next__(self):
        # choose an index
        ix = int(self.prg()) % len(self.l) 
        vx = self.l2[ix]

        # choose a change_balance 
        x1 = self.prg() 
        mdx = self.max_distance_at_index(ix)
        x1 = x1 % mdx 
        s1 = [abs(x1),10**-4]
        x1i = np.argmax(s1) 
        s1_ = [x1,10**-4]
        x1 = s1_[x1i] 

        x2 = -1 if self.c_delta < 0 else 1 
        adj_l = adjust_for_uwpd_change(self.l2,ix,x1*x2,rv=self.srange,\
            d_priority=1,recurse=True) 

        q = self.change_balance
        if type(adj_l) != type(None): 
            self.l2 = adj_l 
            self.change_balance -= (x1 * x2)  
            self.null_ctr = 0
        else: 
            self.null_ctr += 1
        return 

    def max_distance_at_index(self,i):
        d = 0.0
        e0,e1 = 0.0,0.0 
        for j in range(len(self.l)):
            e0 = e0 + abs(self.l[j] - self.srange[0]) 
            e1 = e1 + abs(self.l[j] - self.srange[1]) 
        return max([e0,e1])

    def max_by_index(self):
        return -1 

"""
permuter that changes a base `sequence` according to the wanted 
change in coverage, `coverage_delta` (ratio in [0.,1.]). The 
variables `max_radius` and `super_range` are used to help with 
calculation of coverage. 
"""
class SeqCoveragePermuter(AGV2SeqQualPermuter): 

    def __init__(self,sequence,coverage_delta,max_radius,super_range,prg): 
        super().__init__(sequence,coverage_delta,super_range,prg)

        self.mradius = max_radius 
        self.cov_typeabs,self.cov_typeabs_ = None,None

        self.ocov_typeabs = None 
        self.preproc(self.l)
        return

    def where_is(self,q_value,not_complement:bool=True):

        qx = self.rs if not_complement else self.crs 
        stat0 = is_valid_range(tuple(q_value),False,True)
        stat1 = type(q_value) in {float,np.float32,np.float64}
        assert stat0 or stat1 


        def whereis_cmp(q_):
            if stat1: 
                return q_[0] == q_value or q_[1] == q_value
            return np.all(q_ == q_value)

        for i,q in enumerate(qx):
            if whereis_cmp(q):
                return i 
        return -1 

    def preproc(self,S):
        if type(self.srange) == type(None):
            mn,mx = np.min(S),np.max(S) 
            self.srange = (mn,mx)

        if type(self.mradius) == type(None):
            self.mradius = 0.5  
        assert type(self.mradius) in {float,np.float32,np.float64}

        rs = floatseq_to_rangeseq(S,self.srange,self.mradius)
        self.rs = np.array(to_noncontiguous_ranges(rs))
        self.vci = vector_to_noncontiguous_range_indices(S,self.rs)
        self.crs = np.array(complement_of_noncontiguous_ranges(self.rs,self.srange))
        
        self.update_cov_value() 
        return

    def update_cov_value(self): 
        self.cov_typeabs_ = self.cov_typeabs 
        self.cov_typeabs = range_op(self.rs,default_value=0.,f_inner=np.subtract,f_outer=np.add)

        if type(self.ocov_typeabs) == type(None):
            self.ocov_typeabs = self.cov_typeabs 

    def set_partition(self,m=9): 
        assert type(m) in {int,np.int32,np.int64} 
        assert m > 0 

        # get the new coverage by noncontiguous 
        # ranges
        q = self.srange[1] - self.srange[0]
        q = q * self.c_delta 
        self.new_cov = self.cov_typeabs + q 

        # partition q into m pieces 
        variance = (self.prg() % 10000.) / 10000.
        if variance == 0.0: variance = 0.5 
        num_iter = min([len(self.l) * 4,10**5])
        self.prt = prg_partition_for_float(abs(q),m,self.prg,variance,n=1000,rounding_depth=5)

    def apply_pos_delta(self): 
        for p in self.prt: 
            self.apply_one_pos_delta(p)
            self.format_ranges()

    def apply_one_pos_delta(self,v):

        q = self.pos_delta(v)

        # case: cannot apply delta
        if type(q) in {np.float64,np.float32,float}: 
            return 

        q0_,q1_ = q 
        q0,q1 = q0_ 
        q2,q3 = q1_ 

        # find for `rs`
        ix = self.where_is(np.array(q0),True)
        
        self.rs[ix] = q2 
        ix2 = self.where_is(np.array(q1),False) 
        self.crs[ix2] = q3
        
        self.update_cov_value()
        return

    """
    return:
    - before (range in `rs`,range in `crs`), 
      after (range in `rs`,range in `crs`)
    """
    def pos_delta(self,changean):
        q = self.choose_pos_delta(changean) 

        # case: failure
        if type(q) in {np.float64,np.float32,float}: 
            return q 
            assert False

        n0,n1,ri = q 
        nx = [None,None]
        nx1 = [None,None]
        nxx = [n0,n1] if ri == 0 else [n1,n0] 
        orientation = 0 if nxx[0][1] == nxx[1][0] else 1

        if orientation == 0:
            nx[1] = n0[1] + changean
            nx[0] = n0[0]

            nx1[0] = n0[1] + changean
            nx1[1] = n1[1]
        else:
            nx[0] = n1[0] - changean
            nx[1] = n1[1] 

            nx1[0] = n0[0] 
            nx1[1] = nx[0]

        rx = (n0,n1) if ri == 0 else (n1,n0)
        return rx,(nx,nx1)

    def choose_pos_delta(self,changean):
        if len(self.crs) == 0:
            return 1.0 

        qrs = qualifying_ranges_for_coverage_expansion(self.rs,self.crs,changean,\
            output_type="index")        
        
        # return the maximum 
        if len(qrs) == 0:
            crs_sum = [crs_[1] - crs_[0] for crs_ in self.crs] 
            return max(crs_sum) 

        qi = int(self.prg()) % len(qrs) 
        qiv = qrs[qi]

        # check left and right neighbor 
        neib = neighbors_of_ncrange(self.rs,self.crs,qiv)
        quallr = [False,False] 

        #
        if type(neib[0]) != type(None): 
            dx = neib[0][1] - neib[0][0] 
            if dx >= abs(changean):
                quallr[0] = True 

        if type(neib[1]) != type(None): 
            dx = neib[1][1] - neib[1][0] 
            if dx >= abs(changean):
                quallr[1] = True 

        quallr = np.array([int(q) for q in quallr])
        indices = np.where(quallr == 1)[0] 
        xr = int(self.prg()) % len(indices)
        ix = indices[xr]

        if ix == 0: 
            nr1 = self.rs[qiv] 
            nr0 = neib[0] 
            ref_index=1 
            return nr0,nr1,ref_index 
        else: 
            nr2 = self.rs[qiv] 
            nr3 = neib[1] 
            ref_index= 0 
            return nr2,nr3,ref_index 

    #-------------------------------------------------------------
    
    def change_sequence(self):
        lx = len(self.l)

        seq = []
        for _ in range(lx):
            v2 = self.one_value_in_noncontiguous_range_sequence(False)  
            seq.append(v2) 
        return np.array(seq) 

    def one_value_in_noncontiguous_range_sequence(self,is_complement:bool=False):
        qx = self.crs if is_complement else self.rs 

        lx = [_ for _ in range(len(qx))]
        while len(lx) > 0: 
            i = int(self.prg()) % len(lx) 
            i = lx.pop(i) 

            crange = qx[i] 
            if crange[1] - crange[0] == 0: continue 
            return modulo_in_range(self.prg(),crange)
        
        return None 
    
    def format_ranges(self): 
        self.rs = to_noncontiguous_ranges(self.rs,is_sorted=True)
        self.crs = to_noncontiguous_ranges(self.crs,is_sorted=True)

        def filter_0range(rseq): 
            rs = [] 
            for rs_ in rseq: 
                if rs_[1] - rs_[0] == 0: 
                    continue
                else: 
                    rs.append(rs_) 
            return rs 
        
        self.rs = filter_0range(self.rs) 
        self.crs = filter_0range(self.crs)


    """
    main method 
    """
    def apply(self):
        if self.c_delta > 0: 
            fx = self.apply_pos_delta
        else: 
            fx = self.apply_neg_delta
        fx()

        self.format_ranges() 

        return self.change_sequence()

    def apply_neg_delta(self): 
        for (i,p) in enumerate(self.prt): 
            self.apply_one_neg_delta(p,i) 
            self.format_ranges()

    def apply_one_neg_delta_(self,ncrange,ncrange_comp,v):
        # get orientation
        o = 0 
        if ncrange[0] == ncrange_comp[1]:
            o = 1 
        else: 
            assert ncrange[1] == ncrange_comp[0] 
        
        if o == 0:
            ncrange[1] -= abs(v)
            ncrange_comp[0] = ncrange[1] 
        else: 
            ncrange[0] += abs(v)
            ncrange_comp[1] = ncrange[0]  
        return ncrange,ncrange_comp 


    def default_random_shrink_(self,dx=2.0):
        ri = int(self.prg()) % len(self.rs) 
        rx = deepcopy(self.rs[ri] )
        d = (rx[1] - rx[0]) / dx
        rx2 = deepcopy(rx)
        rx2[1] = rx2[0] + d 

        d2 = rx[1] - rx2[1] 
        neib = neighbors_of_ncrange(self.rs,self.crs,ri)[1] 
        self.rs[ri] = rx2

        if type(neib) != type(None):
            ix = self.where_is(tuple(neib),False) 
            self.crs[ix] = np.array([rx2[1],self.crs[ix][1]]) 
        return ri,d2 

    def apply_one_neg_delta(self,v,i):
        ri = qualifying_ranges_for_coverage_shrinkage(self.rs,v,output_type="index")

        # no qualifying ranges for partition. have to distribute
        # error to the rest
        if len(ri) == 0: 
            ri,ex = self.default_random_shrink_() 
            ox = self.prt[i] - ex 
            self.prt[(i+1)%len(self.prt)] += ex  
            self.update_cov_value()
            return 
        assert len(ri) > 0 

        rii = int(self.prg()) % len(ri)
        xr = ri[rii] 

        # get complementary neighbors 
        neib = neighbors_of_ncrange(self.rs,self.crs,xr)

        possible_complement_ranges = []

        if type(neib[0]) != type(None): 
            if abs(v) <= neib[0][1] - neib[0][0]:
                possible_complement_ranges.append(neib[0]) 

        if type(neib[1]) != type(None): 
            if abs(v) <= neib[1][1] - neib[1][0]:
                possible_complement_ranges.append(neib[1]) 

        rii2 = int(self.prg()) % len(possible_complement_ranges)
        crange = possible_complement_ranges[rii2] 
        ix = self.where_is(crange,False)

        r0 = self.rs[xr] 
        r1 = self.crs[ix] 

        q0,q1 = self.apply_one_neg_delta_(r0,r1,v)
        self.rs[xr] = q0
        self.crs[ix] = q1 
        self.update_cov_value()
