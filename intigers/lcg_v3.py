from .lgc_v2 import * 
from .tvec import * 

class LCGV3(LCGV2): 

    def __init__(self,start,m,a,n0,n1,sc_size:int,preproc_gd:bool=False):
        super().__init__(start,m,a,n0,n1,sc_size,preproc_gd)
        self.tv = None
        self.

    def __next__(self):
        return -1

    def set_tv(self,tv):
        assert type(tv) == TrinaryVec
        self.tv = tv 
        return

    def autoset_tv(self,gen_type,l,ext_prg=None):
        assert gen_type in {1,2}

        qx = ext_prg if type(ext_prg) != \
            None else self.__next__ 

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