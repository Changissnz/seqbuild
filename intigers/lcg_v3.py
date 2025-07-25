from .lcg_v2 import * 
from .tvec import * 
from morebs2.numerical_generator import modulo_in_range

"""
"""
def is_value_below_modulo_range_length(v,mr):
    assert is_valid_range(mr,True,False)  
    return abs(v) <= mr[1] - mr[0]

# NOTE: not the most efficient approach
def iterative_adjustment_for_direction__modrange(pv,av,sign,modrange,\
    modindex,moddir):  
    assert is_valid_range(modrange,True,False) 

    x1 = modulo_in_range(av,modrange)
    
    # case: do nothing
    if to_trinary_relation(x1,pv) == sign: 
        return modrange

    newrange = [modrange[0],modrange[1]]

    stat = True
    while stat:
        newrange[modindex] += moddir 
        x1 = modulo_in_range(av,newrange)
        stat = to_trinary_relation(x1,pv) != sign
    return newrange 

"""
NOTE: Method is specifically for `av` that does not satisfy the 
`is_value_below_modulo_range_length(av,modrange)` relation. 

Adjusts a modulo range `modrange` to M_r so that the output from 
`modulo_in_range(av,M_r)` satisfies the `sign`, -1 for `pv > av` or 
1 for `pv < av`. 

pv := "parent" value, also a "reference" for the trinary relation with 
      `av`.
av := actual value, v = m * pv + a. 
sign := wanted binary relation between `av` and `pv`. 
modrange := base modulo range to adjust. 
"""
def multimodular_number__modulo_range_adjustment(pv,av,sign,modrange):
    assert sign in {-1,1}

    if sign == 1: 
        md,mi = 1,1
    else: 
        md,mi = -1,0
    return iterative_adjustment_for_direction__modrange(pv,av,sign,\
        modrange,mi,md) 

"""
NOTE: Method is specifically for `av` that satisfies the 
`is_value_below_modulo_range_length(av,modrange)` relation. 

wv := wanted value
av := actual value
cv := current value, the output from `modulo_in_range(av,modrange)`
modrange := base modulo range to adjust. 
"""
def unimodular_number__modulo_range_adjustment(wv,av,cv,modrange):

    sign = to_trinary_relation(wv,cv) 
    assert sign != 0 

    # adjust modrange if necessary
    cv_ = cv 
    mr = [modrange[0],modrange[1]] 
    if sign == 1: 
        if mr[1] < wv:
            mr[1] = wv
            cv_ = modulo_in_range(av,mr)
            mr[1] += wv - cv_ 
        else: 
            if mr[1] == wv: 
                mr[1] = wv + 1 
            d = wv - av 
            mr[0] = d 
    else:
        if mr[0] > wv: 
            mr[0] = wv
            cv_ = modulo_in_range(av,mr)
            mr[1] = mr[1] - cv + wv 
        elif mr[0] == wv: 
            mr[1] = mr[1] - (cv_ - mr[0])
        else: 
            mr[0] -= (cv_ - wv) 
    return mr 

"""
Calculates a new modrange M_r such that for the actual 
value `av`, `modulo_in_range(av,M_r) == pv`. 
"""
def modrange_for_congruence(pv,av,modrange):
    assert type(pv) in {int,np.int32,np.int64}
    assert type(av) in {int,np.int32,np.int64}
    assert is_valid_range(modrange,True,False) 

    if not is_value_below_modulo_range_length(av,modrange): 
        if pv >= 0: 
            mr2 = [0,abs(av)] 

            m = 1 if av < 0 else -1 
            mr2[1] += m * pv 
        else:            
            av_ = -av if av > 0 else av 

            pvs = pv >= 0 
            avs = av >= 0 
            s = -1 if pvs == avs else 1 
            mr2 = [av_ + 1,0] 
            mr2[0] = mr2[0] + s * pv - 1
        return mr2 

    cv = modulo_in_range(av,modrange)
    difference = cv - pv 
    mr2 = [modrange[0],modrange[1]]
    
    if av >= 0: 
        mr2[0] -= difference 
    else: 
        mr2[1] -= difference 
    return mr2 

#--------------------------------------------------------------------- 

class LCGV3(LCGV2): 

    def __init__(self,start,m,a,n0,n1,sc_size:int,preproc_gd:bool=False,\
        prg=None):
        super().__init__(start,m,a,n0,n1,sc_size,preproc_gd)
        self.tv = None
        self.tvi = None 
        self.r_ = [self.r[0],self.r[1]]
        self.prg = prg

    def __next__(self):
        did_fire = self.fired
        s_ = self.s 
        q = super().__next__()

        # case: first element, no trinary vec comparison 
        if not did_fire:
            return q

        ## print("S_: ",s_, " Q: ",q, " TVI: ",self.tvi, " SG: ",self.tv[self.tvi])

        if type(self.tv) == TrinaryVec: 
            assert type(self.prg) != type(None)

            sc = self.tv[self.tvi]
            mr = self.adjust_modulo_range(s_,q,sc,self.prg)
            if type(mr) != type(None): 
                self.r = [int(mr[0]),int(mr[1])]
                q = s_ * self.m + self.a
                q = modulo_in_range(q,self.r) 
            self.tvi = (self.tvi + 1) % len(self.tv) 
            self.s = q 
        return q 

    def adjust_modulo_range(self,reference,new_value,sign_change,prg2):
        # case: sign already satisfied
        if to_trinary_relation(new_value,reference) == sign_change: 
            return None 

        actual_value = self.m * reference + self.a 

        # case: congruence 
        if sign_change == 0: 
            r_ = modrange_for_congruence(reference,actual_value,self.r)
            return r_ 

        stat = is_value_below_modulo_range_length(actual_value,self.r) 
        if stat: 
            wv = reference 
            d = modulo_in_range(prg2(),[0,self.r[1] - self.r[0]])
            wv += (d * sign_change)
            return unimodular_number__modulo_range_adjustment(wv,actual_value,\
                new_value,self.r)

        return multimodular_number__modulo_range_adjustment(reference,\
            actual_value,sign_change,self.r)

    def set_tv(self,tv):
        assert type(tv) == TrinaryVec
        self.tv = tv 
        self.tvi = 0 
        return

    def autoset_tv(self,gen_type,l,ext_prg=None):
        assert gen_type in {1,2}

        qx = ext_prg if type(ext_prg) != \
            type(None) else self.__next__ 

        tv = None
        if gen_type == 1:
            L = [-1,0,1]
            bv = L[qx() % 3]
            cr = (qx() % 1000) / 1000.0 
            tv = TrinaryVec.one_instance__v1(bv,l,cr,qx) 
        else: 
            k0 = (qx() % 1000) / 1000.0
            k1 = (qx() % 1000) / 1000.0
            k2 = (qx() % 1000) / 1000.0
            s = k0 + k1 + k2 
            k0,k1,k2 = zero_div0(k0,s),\
                zero_div0(k1,s),zero_div0(k2,s)
            fm = {-1:k0,0:k1,1:k2}
            tv = TrinaryVec.one_instance__v2(l,fm,qx)
        self.set_tv(tv)

    def set_prg(self,prg):
        assert type(prg) in {FunctionType,MethodType}
        self.prg = prg 