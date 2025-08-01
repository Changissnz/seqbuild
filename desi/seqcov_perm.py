from mini_dm.ag_ext import * 
from morebs2.numerical_generator import prg_partition_for_float

class SequenceCoveragePermuter: 

    def __init__(self,sequence,coverage_delta,max_radius,super_range,prg): 
        assert is_vector(sequence) 
        assert type(coverage_delta) in {float,np.float32,np.float64}
        assert -1. <= coverage_delta <= 1. 
        #assert type(max_radius) in {float,np.float32,np.float64}
        assert is_valid_range(super_range,False,False) or type(super_range) == type(None)
        assert type(prg) in {MethodType,FunctionType} 
        
        self.l = sequence 
        self.coverage_delta = coverage_delta 
        self.mradius = max_radius 
        self.srange = super_range 
        self.prg = prg  

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
        q = q * self.coverage_delta 
        self.new_cov = self.cov_typeabs + q 

        # partition q into m pieces 
        variance = (self.prg() % 10000.) / 10000.
        if variance == 0.0: variance = 0.5 
        num_iter = min([len(self.l) * 4,10**5])
        self.prt = prg_partition_for_float(abs(q),m,self.prg,variance,n=1000,rounding_depth=5)

    def apply_pos_delta(self): 
        for p in self.prt: 
            self.apply_one_pos_delta(p)

    def apply_one_pos_delta(self,v):
        q0_,q1_ = self.pos_delta(v)
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
        if type(q) == float: 
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

    def apply_neg_delta(self): 
        for (i,p) in enumerate(self.prt): 
            self.apply_one_neg_delta(p,i) 


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
            ix = self.where_is(tuple(neib)) 
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
