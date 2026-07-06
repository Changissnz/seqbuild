from .seq_struct import * 
from morebs2.numerical_generator import prg_choose_n,modulo_in_range
from morebs2.matrix_methods import is_valid_range,vector_to_string,subvec 
from morebs2.measures import to_trinary_relation
from .extraneous import round_to_trinary_vector,prg__single_to_int,prg__single_to_trinary_vector
from types import MethodType,FunctionType


class TrinaryVec(IntSeq):

    def __init__(self,l):
        assert set(l).issubset({0,1,-1})
        super().__init__(l)
        self.index = 0 

    def __next__(self):
        q = super().__next__() 
        self.index = self.index % len(self.l) 
        return q 
    
    @staticmethod
    def omit_value(tv,target_value,prg):
        assert type(tv) == TrinaryVec 
        assert target_value in {-1,0,1} 
        assert type(prg) in {MethodType,FunctionType} 

        lx = []

        repset = {-1,0,1} - {target_value}
        repset = sorted(repset) 
        for l_ in tv.l:
            if l_ == target_value: 
                i = int(prg()) % 2 
                r = repset[i]
                lx.append(r)
            else: 
                lx.append(l_) 
        return TrinaryVec(lx) 

    # TODO: test these two methods
    @staticmethod
    def one_instance__v1(base_value,length,change_ratio,prg): 
        assert base_value in {-1,0,1}
        assert change_ratio >= 0.0 and change_ratio <= 1.0 
        
        prg_ = prg__single_to_int(prg) 
        vec = np.zeros((length,),dtype=np.int32) \
            + base_value 

        n = int(ceil(change_ratio * length))
        Q = [i for i in range(length)]
        ix = prg_choose_n(Q,n,prg_,is_unique_picker=True)

        L = {-1,0,1} - set([base_value]) 
        L = sorted(L) 

        for i in ix:
            d = L[prg_() % 2]
            vec[i] = d 
        return TrinaryVec(vec)

    @staticmethod
    def one_instance__v2(length,frequency_map,prg): 
        s = sum(frequency_map.values())
        assert abs(1.0 - s) < 10 ** -5, "got {}".format(frequency_map)
        assert set(frequency_map.keys()).issubset({-1,0,1})
        rx = [(k,int(round(v * length))) for k,v in frequency_map.items()]
        rx = np.array(rx,dtype=np.int32) 
        sx = length - sum(rx[:,1])

        # choose random index to add calibration to 
        if sx != 0:
            rp = int(prg()) % len(rx) 
            rx[rp,1] = rx[rp,1] + sx 

        Q = [i for i in range(length)]
        vec = np.zeros((length,),dtype=np.int32)
        prg_ = prg__single_to_int(prg)
        for q in rx:
            k,v = q[0],q[1] 
            if v == 0: 
                ##print("...strange...")
                break 

            q2 = prg_choose_n(Q,v,prg_,is_unique_picker=True)
            vec[q2] += k
        return TrinaryVec(vec)

    @staticmethod
    def one_instance__v3(l,prg,super_range): 
        assert is_valid_range(super_range,False,False) 
        
        prg_ = prg__single_to_int(prg)

        lx = [modulo_in_range(prg_(),super_range) for \
            _ in range(l)]
        V = round_to_trinary_vector(lx,is_distance_roundtype=False)
        return TrinaryVec(V) 

    """
    calculates rightward contiguous range for l[i]. 
    """
    def contiguous_range(self,i):
        q = self.l[i]
        start,end = i,None
        for j in range(i+1,len(self.l)): 
            if self.l[j] != q: 
                end = j
                break 

        if type(end) == type(None):
            end = len(self.l) 
        return [start,end] 


#--------------------------------------------------

def std_trinary_vec__oscillating(ndim): 
    s = ceil(ndim / 3)
    q = [-1,0,1] * s
    return np.array(q[:ndim]) 

def std_trinary_vec__two_halves(ndim,half1,half2): 
    assert half1 != half2 
    assert half1 in {-1,0,1}
    assert half2 in {-1,0,1} 

    l = ceil(ndim / 2) 
    q0 = [half1] * l 
    q1 = [half2] * (ndim - l) 
    return np.array(q0 + q1) 

"""
generates 8 trinary vectors derived from `V`, 
a trinary vector. These trinary vectors are 
for the goal of (u)nique (o)bjective.  
"""
def trinary_vec_delta_sequence__type_UO(V):
    assert set(V).issubset({-1,0,1}) 

    #S0_ = np.ones((len(V),))
    S0 = std_trinary_vec__oscillating(len(V))
    S1 = std_trinary_vec__two_halves(len(V),-1,0)
    S2 = std_trinary_vec__two_halves(len(V),0,1)

    #yield modulo_in_range(V + S0_,[-1,2])
    #yield modulo_in_range(V * -S0_,[-1,2])

    yield modulo_in_range(V + S0,[-1,2])  
    yield modulo_in_range(V + S1,[-1,2]) 

    yield modulo_in_range(V - S0,[-1,2])  
    yield modulo_in_range(V + S2,[-1,2])  

    l = ceil(len(V) / 3)
    
    yield modulo_in_range(V + subvec(S0,l,len(V)),[-1,2])  
    yield modulo_in_range(V + subvec(S1,l,len(V)),[-1,2]) 

    yield modulo_in_range(V - subvec(S0,l,len(V)),[-1,2])  
    yield modulo_in_range(V + subvec(S2,l,len(V)),[-1,2])  

# NOTE: this is an attempt-to-generate process; 
#       not guaranteed to produce the wanted `m`
#       unique trinary vectors. For the `m` wanted 
#       vectors, makes ceil(`attempt_ratio`` * `m`) 
#       attempts. 
# NOTE: do NOT call this function for large parameters 
#       `ndim` and `m`, since all generated trinary vectors 
#       are stored in memory.  
def generate_m_unique_trinary_vectors(ndim,m,prg,attempt_ratio=2.0):  

    if m > 3 ** ndim: 
        m = 3 ** ndim
        
    num_attempts = ceil(attempt_ratio * m)

    L = set() 
    prev_unique = None 
    prg_ = prg__single_to_trinary_vector(prg,ndim)

    while num_attempts > 0 and m > 0: 
        l = prg_() 
        x = vector_to_string(l,int) 
        
        # case: unique
        if x not in L:
            L |= {x}
            prev_unique = l 
            num_attempts -= 1 
            m -= 1 
            yield l
        else: 
            # try looking for a unique one through derivation sequence on 
            # l 
            Q = trinary_vec_delta_sequence__type_UO(prev_unique) 
            for q in Q: 
                x = vector_to_string(q,int) 
                num_attempts -= 1
                if x not in L: 
                    L |= {x}
                    prev_unique = q 
                    m -= 1  
                    yield q 
                    break 
            
# standalone function, outside of class<LCGV3> that also employs an adjustment scheme 
# for this same objective: adjusting the values of v, if need be, to suit the trinary vector t 
def ternary_adjustment(v,t,prg):
    assert len(v) - 1 == len(t)

    if len(v) < 2: return v 

    def next_from_prg(v1_,r): 
        x = abs(prg())

        if x == 0: x = 1 
        return x * r 

    def adjust(v0_,v1_,r): 
        c = 50 
        while to_trinary_relation(v1_,v0_) != r: 
            q = next_from_prg(v1_,r) 
            v1_ = v1_ + q            

        return v1_ 
    
    v_ = np.copy(v) 

    for i in range(1,len(v_)):
        v0,v1 = v_[i-1],v_[i]  
        R = to_trinary_relation(v1,v0)
        R_ = t[i - 1] 
        if R == R_: continue 

        if R_ == 0: 
            v_[i] = v0 
            continue 

        v2 = adjust(v0,v1,R_) 
        v_[i] = v2 
    return v_  





