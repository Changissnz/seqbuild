"""
Planned file that works on vector-singleton fits between vectors. 
"""

from morebs2.matrix_methods import is_number,is_vector,is_valid_range
from morebs2.numerical_generator import modulo_in_range
from morebs2.measures import zero_div 
from intigers.extraneous import to_trinary_relation_v2,safe_div
import numpy as np 

DEFAULT_AFFINEVEC_DIMRANGE = [3,17]

#-------------------- some methods on ratios 

# parameters for method<ratio_vector> 
ratio_vector__types = ["a","s"]
ratio_vector__parameter = {"a": ["min 1.0","max 1.0"],"s":[0,1]} 
    # active var only for "a" 
ratio_vector__parameter2 = [-1,0,1]

"""
outputs the non-1.0 value from the scaling of 
pair (q0,q1) according to pair relation of 
min|max 1.0. 
"""
def ratio__type_asymmetric(q0,q1,vec_type):
    assert vec_type in {"min 1.0", "max 1.0"}

    if vec_type == "min 1.0":
        m0,m1 = min([q0,q1]),max([q0,q1])
    else: 
        m0,m1 = max([q0,q1]),min([q0,q1])
    return zero_div(m1,m0,0.0)

"""
outputs the corresponding pair (R(q0),R(q1)); 
The set {R(q0),R(q1)} has the element 1.0 and either 
a < 1.0 or > 1.0 value. The resulting output's index at 
1.0 makes that index the reference index.

NOTE: method is the complete output from the calculation 
in method<ratio__type_asymmetric>.
"""
def ratio__type_asymmetric__v2(q0,q1,vec_type):
    q = ratio__type_asymmetric(q0,q1,vec_type)

    x = [None,None]
    if vec_type == "min 1.0": 
        i0 = np.argmin([q0,q1])
    else: 
        i0 = np.argmax([q0,q1])
    i1 = (i0 + 1) % 2 
    x[i0] = 1.0
    x[i1] = q 

    return x 

def ratio__type_symmetric(q0,q1,ref=0):
    assert ref in {0,1}, "got {}".format(ref)

    x = abs(q0) + abs(q1) 

    if ref == 0: 
        return zero_div(q0,x,0.0)
    return zero_div(q1,x,0.0)

"""
Calculates a ratio vector `q0/q1` that is either a float 
or an n-dimensional vector. The parameter space is the following: 

q0 := float|vector,numerator 
q1 := float|vector,denumerator 
rtype := (s)ymmetric or (a)symmetric 
parameter := primary parameter used; 
            if `s` -> [0,1], 
            if `a` -> [min 1.0,max 1.0].
parameter2 := used if rtype is `s`; 
                [-1,0,1]
parameter3 := dict, used for `auto` rtype. 
                key is `s`|`a`,
                value is default argument. 
auto_output := 0 for ratio vector, 
                1 for ratio vector + corresponding "s/a" status 
"""
def ratio_vector(q0,q1,rtype,parameter,parameter2,parameter3=None,\
    auto_output=0):
    assert rtype in {"a","s","auto"}
    if rtype in {"a","auto"}:
        assert parameter2 in {-1,0,1}

    if rtype == "auto":
        assert set(parameter3.keys()) == set(["a","s"])
        assert parameter3["a"] in ratio_vector__parameter["a"]
        assert parameter3["s"] in ratio_vector__parameter["s"]

    def qf__autolabel(qx0,qx1,param):
        # case: use asymmetric labeling
        if (qx0 < 0.0 and qx1 > 0.0) or \
            (qx0 > 0.0 and qx1 < 0.0): 
            x = ratio__type_asymmetric(qx0,qx1,parameter3["a"])
            if auto_output == 0: 
                return x 
            return x,"a" 
        x = ratio__type_symmetric(qx0,qx1,parameter3["s"])
        return x if auto_output == 0 else (x,"s")
    
    if rtype == "auto":
        qf = qf__autolabel
    elif rtype == "a":
        def qf(qx0,qx1,param): 
            if parameter2 == -1: 
                return ratio__type_asymmetric(qx0,qx1,param)
            xr = ratio__type_asymmetric__v2(qx0,qx1,param) 
            return xr[parameter2] 
    else: 
        qf = ratio__type_symmetric

    i = 0 if is_vector(q0) else None
    j = 0 if is_vector(q1) else None 

    if i != 0 and j != 0:
       return qf(q0,q1,parameter)

    if i == j and i == 0:
        assert len(q0) == len(q1) 

    lx = [] 
    stat1,stat2 = True,True 
    while stat1 and stat2: 

        x1 = q0[i] if type(i) == int else q0 
        x2 = q1[j] if type(j) == int else q1 

        lx_ = qf(x1,x2,parameter)
        lx.append(lx_)

        stat1,stat2 = True,True

        if type(i) == int: 
            i += 1 
            stat1 = not (i >= len(q0))
        if type(j) == int: 
            j += 1 
            stat2 = not (j >= len(q1))
    
    return np.array(lx) 

#--------------------------- 

DEFAULT_MA_DISTANCE_FUNCTION = lambda x,x2: np.abs(x) + np.abs(x2)

"""
container that describes the geometric relation of 
an <AffineDelta>'s variables.
"""
class MADescriptor:

    def __init__(self,rv_vec,rvt_vec,t_vec,s_vec,d_vec,ma_order=None):
        # ratio vector 
        self.rv_vec = rv_vec 
        # symmetric status vector 
        self.rvt_vec = rvt_vec
        # relation of first w/ second (according to MA order)
        self.t_vec = t_vec
        # sign vector 
        self.s_vec = s_vec
        # absolute distance vector 
        self.d_vec = d_vec 
        self.ma_order = ma_order
        return 
    
    def __str__(self):
        s = ""
        s += "* RV: " + str(self.rv_vec) + "\n"
        s += "* RVT: " + str(self.rvt_vec) + "\n"
        s += "* T: " + str(self.t_vec) + "\n"
        s += "* S: " + str(self.s_vec) + "\n"
        s += "* D: " + str(self.d_vec) + "\n"
        return s
    
    def dim(self):
        if not is_vector(self.rv_vec): 
            l = 0  
        else: 
            l = len(self.rv_vec)
        return l 

    # outputs (multiple m,adder a) pair; 
    # m,a are each one of float or vector. 
    def solve(self,ma_dim): 
        d = self.dim()

        assert max(ma_dim) == d

        X = [None,None]
         
        X[0] = 0.0 if ma_dim[0] == 0 else \
            np.zeros((ma_dim[0],),dtype=float)

        X[1] = 0.0 if ma_dim[1] == 0 else \
            np.zeros((ma_dim[1],),dtype=float)
        
        def add_at_i(i,y0,y1):

            if is_number(X[0],set()): 
                X[0] = y0
            else: 
                X[0][i] = y0 

            if is_number(X[1],set()): 
                X[1] = y1
            else: 
                X[1][i] = y1
            return        

        # single for both 
        if d == 0: 
            return self.solve_at_(\
                self.rvt_vec,self.d_vec,\
                self.s_vec,self.rv_vec,self.t_vec)

        for i in range(d): 
            q0,q1 = self.solve_at(i)
            add_at_i(i,q0,q1)
        return X 
    
    def solve_at(self,i):
        def varvalue(value,i):
            if is_number(value,set()): return value 
            return value[i] 

        a_stat = varvalue(self.rvt_vec,i) 
        d = varvalue(self.d_vec,i)
        s = varvalue(self.s_vec,i)
        r = varvalue(self.rv_vec,i)
        t = varvalue(self.t_vec,i)
        return self.solve_at_(a_stat,d,s,r,t) 
    
    def solve_at_(self,a_stat,d,s,r,t):
    
        x0,x1 = None,None 
        if a_stat == 'a': 
            # min 1.0
            ##a_stat = "min 1.0" if q <= 1.0 else "max 1.0"
            total = 1.0 + r
            base_one = safe_div(1.0,total)
            other = 1.0 - base_one

            base_share = d * base_one
            other_share = d * other 

            qshare = sorted([base_share,other_share])
            index = 1 if t == 1 else 0 

            if s == -1:
                qshare[index] = -1 * qshare[index] 
            else: 
                qshare[(index+1)%2] = -1 * \
                    qshare[(index+1)%2] 

            q0 = qshare[index] 
            q1 = qshare[(index + 1) % 2]
            return (q0,q1)
        
        q = r * d 
        q2 = d - q 
        
        if s == 0: s = 1 

        q,q2 = q * s,\
            q2 * s
        return q,q2 
    
    @staticmethod
    def from_AffineDelta(ad,p3,d_operator=lambda x2,x1:x2-x1):

        q0,q1 = ad.m,ad.a
        if ad.ma_order == 0:
            qx0 = [ad.m,ad.a]
        else: qx0 = [ad.a,ad.m]

        parameter2 = -1
        rv = ratio_vector(qx0[0],qx0[1],"auto",None,\
            parameter2,parameter3=p3,auto_output=1)
        if ad.ma_order == 0:
            qref,qref2 = ad.m,ad.a 
        else: 
            qref,qref2 = ad.a,ad.m
        rx = [qref,qref2]

        tvec = to_trinary_relation_v2(np.abs(rx[0]),np.abs(rx[1]))
        qref3 = 0.0 if not is_vector(rx[0]) else np.zeros((len(rx[0]),))
        svec = to_trinary_relation_v2(qref,qref3) 

        dvec = np.abs(d_operator(qref2,qref)) 
        
        if type(rv) == tuple: 
            qx1,qx2 = np.abs(rv[0]),rv[1]
        else: 
            qx1,qx2 =np.abs(np.array(rv[:,0],dtype=float)),\
                rv[:,1]
        return MADescriptor(qx1,qx2,tvec,\
            svec,dvec,ad.ma_order) 

"""
designed for only vector and singleton values. 
This is the primary structure that fits between

X_{input} <-> X_{output}; 
X_{input},X_{output} each of type vector or singleton. 
"""
class AffineDelta: 

    def __init__(self,m,a,ma_order:int): 
        assert is_number(m,set()) or is_vector(m) 
        assert is_number(a,set()) or is_vector(a) 
        assert ma_order in {0,1}

        self.m = m
        self.a = a 
        self.ma_order = ma_order
        return
    
    def type(self):
        q = int(is_vector(self.m))
        q2 = int(is_vector(self.a)) 
        return q,q2 

    
    def __str__(self):
        s = "m: {}".format(self.m)
        s2 = "a: {}".format(self.a)
        o = "o: {}".format(self.ma_order)
        return s + "\n" + s2 + "\n" + o + "\n"

    #------------------------------------- i/o methods 

    def fit(self,x):
        if self.ma_order == 0:
            return x * self.m + self.a
        return (x + self.a) * self.m 
    
    def op1(self,x):
        if self.ma_order: return x * self.m 
        return x + self.a 

    def diff(self,x,x2): 
        return self.fit(x) - self.fit(x2) 
    
    def expected_diff(self,target_value,x,op_type):
        assert op_type in {"one","all"}

        if op_type == "one": 
            return target_value - self.op1(x)
        return target_value - self.fit(x) 
    
    def to_ma_descriptor(self,p3,d_operator=lambda x,x2:x2 - x):
        return MADescriptor.from_AffineDelta(self,p3,d_operator)

    def delta(self,dfunc): 
        m,a = dfunc(self.m,self.a) 
        return AffineDelta(m,a,self.ma_order)
    
    '''
    contribution vector for M and A w.r.t. input x 
    '''
    def cvec(self,x):

        q1 = self.op1(x)
        q2 = self.fit(x)  

        t1 = np.abs(x - q1)
        t2 = np.abs(q2 - q1)
        t = t1 + t2  

        rx = safe_div(t1,t)
        rx2 = safe_div(t2,t)
        return rx,rx2 

    #------------------------------ instantiation methods      

    @staticmethod
    def one_instance_(prg,r_out1,r_out2,dim_range=None,ma_order=None):
        dtypes = ["vec","float"]

        m_type = dtypes[prg() % 2]
        a_type = dtypes[prg() % 2]
        return AffineDelta.one_instance(m_type,a_type,prg,r_out1,r_out2,\
            dim_range,ma_order)

    @staticmethod 
    def one_instance(m_type,a_type,prg,r_out1,r_out2,dim_range=None,ma_order=None):
        M,A = AffineDelta.one_instance_vars(m_type,a_type,prg,r_out1,r_out2,dim_range)

        if type(ma_order) != type(None):
            assert ma_order in {0,1}
        else: 
            ma_order = prg() % 2 

        return AffineDelta(M,A,ma_order)        
    
    @staticmethod 
    def one_instance_vars(m_type,a_type,prg,r_out1,r_out2,dim_range=None):
        assert m_type in {"vec","float"}
        assert a_type in {"vec","float"}

        if (m_type == "vec" or a_type == "vec") and \
            type(dim_range) == type(None): 
            dim_range = DEFAULT_AFFINEVEC_DIMRANGE
        #else: 
        #    assert is_valid_range(dim_range)

        def output_one_value(value_idn,sz=None):
            assert value_idn in {0,1} 

            vtype = None
            if value_idn == 0:
                rx = r_out1() 
                vtype = m_type
            else: 
                rx = r_out2()
                vtype = a_type
                    
            if vtype == "vec":
                if type(sz) == type(None): 
                    ql = modulo_in_range(prg(),dim_range) 
                else: 
                    ql = sz 

                q = np.zeros((ql,),dtype=np.float64) 
                assert is_valid_range(rx,False)

                for i in range(len(q)):
                    v = modulo_in_range(prg(),rx) 
                    q[i] = v 
            else:
                q = modulo_in_range(prg(),rx) 
            return q 


        M = output_one_value(0)
        sz = None
        if m_type == "vec":
            sz = len(M) 
        A = output_one_value(1,sz)  

        return M,A 