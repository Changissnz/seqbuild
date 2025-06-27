"""
Planned file that works on vector-singleton fits between vectors. 
"""

from morebs2.matrix_methods import is_number,is_vector,is_valid_range
from morebs2.numerical_generator import modulo_in_range
from morebs2.measures import zero_div 
import numpy as np 

DEFAULT_AFFINEVEC_DIMRANGE = [3,17]


#-------------------- some methods on ratios 

# parameters for method<ratio_vector> 
ratio_vector__types = ["a","s"]
ratio_vector__parameter = {"a": ["min 1.0","max 1.0"],"s":[0,1]} 
    # active var only for "a" 
ratio_vector__parameter2 = [-1,0,1]

def ratio__type_asymmetric(q0,q1,vec_type):
    assert vec_type in {"min 1.0", "max 1.0"}

    if vec_type == "min 1.0":
        m0,m1 = min([q0,q1]),max([q0,q1])
    else: 
        m0,m1 = max([q0,q1]),min([q0,q1])
    return zero_div(m1,m0,0.0)

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

# parameter3 := used for auto (default for a, default for s)
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
       print("-- output")
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
    
    def ma_ratio(self,ratio_type,parameter):
        assert ratio_type in {"a","s"}

        """
        rv = ratio_vector(self.m,self.a,\
            ratio_type,parameter,parameter2):
    assert rtype in {"a","s"}

    if rtype == "a":
        assert parameter2 in {-1,0,1}


        if ratio_type == "a": 
            return 

        ###### 
def safe_div(V1,V2):
def ratio__type_asymmetric(q0,q1,vec_type):
    assert vec_type in {"min 1.0", "max 1.0"}

    if vec_type == "min 1.0":
        m0,m1 = min([q0,q1]),max([q0,q1])
    else: 
        m0,m1 = max([q0,q1]),min([q0,q1])
    return zero_div(m1,m0,0.0)

def ratio__type_symmetric(q0,q1,ref=0):
        return -1 
        """ 

    def delta(self,dfunc): 
        m,a = dfunc(self.m,self.a) 
        return AffineDelta(m,a,self.ma_order)
    


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