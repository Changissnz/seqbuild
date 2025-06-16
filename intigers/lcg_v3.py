from .lcg_v2 import * 
from .tvec import * 
from morebs2.numerical_generator import modulo_in_range

def is_value_below_modulo_range_length(v,mr):
    assert is_valid_range(mr,True,False)  
    return abs(v) <= mr[1] - mr[0]

# NOTE: not the most efficient approach
def iterative_adjustment_for_direction__modrange(pv,av,sign,modrange,\
    modindex,moddir):  
    assert is_valid_range(modrange,True,False) 

    x1 = modulo_in_range(av,modrange)
    
    # case: do nothing
    if to_trinary_relation(pv,x1) == sign: 
        return modrange

    newrange = [modrange[0],modrange[1]]
    ##q = -1 if modindex == 0 else 1
    #moddir = None 
    #if modindex == 0 and av > 

    stat = True
    while stat:
        newrange[modindex] += moddir 
        x1 = modulo_in_range(av,newrange)
        stat = to_trinary_relation(pv,x1) != sign
    return newrange 

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
            if av > 0: 
                av_ = -av 
            else: 
                av_ = av 

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

class LCGV3(LCGV2): 

    def __init__(self,start,m,a,n0,n1,sc_size:int,preproc_gd:bool=False):
        super().__init__(start,m,a,n0,n1,sc_size,preproc_gd)
        self.tv = None
        self.r_ = [self.r[0],self.r[1]]

    def __next__(self):
        return -1

    def set_tv(self,tv):
        assert type(tv) == TrinaryVec
        self.tv = tv 
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
        self.tv = tv 