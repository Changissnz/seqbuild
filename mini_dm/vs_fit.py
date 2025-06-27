"""
Planned file that works on vector-singleton fits between vectors. 
"""

from morebs2.matrix_methods import is_number,is_vector,is_valid_range
from morebs2.numerical_generator import modulo_in_range
from morebs2.measures import zero_div 
import numpy as np 

DEFAULT_AFFINEVEC_DIMRANGE = [3,17]


"""
most under or least over; 

takes a delta function to increment point in search 
"""
def mulo___inc1(p_current,p_ref,diff_func,delta_func,term_func,max_iterations:int):
    d = diff_func(p_ref,p_current)
    q = delta_func(p_current)
    stat = term_func(q) 

    return d,q,stat

def ratio__type_asymmetric(q0,q1,vec_type):
    assert vec_type in {"min 1.0", "max 1.0"}

    if vec_type == "min 1.0":
        m0,m1 = min([q0,q1]),max([q0,q1])
    else: 
        m0,m1 = max([q0,q1]),min([q0,q1])
    return zero_div(m1,m0,0.0)

def ratio__type_symmetric(q0,q1,ref=0):
    assert ref in {0,1}

    x = abs(q0) + abs(q1) 

    if ref == 0: 
        return zero_div(q0,x,0.0)
    return zero_div(q1,x,0.0)

"""
designed for only vector and singleton values 
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
    
    def diff(self,x,x2): 
        return self.fit(x) - self.fit(x2) 
    
    def delta(self,dfunc): 
        m,a = dfunc(self.m,self.a) 
        return AffineDelta(m,a,self.ma_order)
    
    def op1(self,x):
        if self.ma_order: return x * self.m 
        return x + self.a 

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