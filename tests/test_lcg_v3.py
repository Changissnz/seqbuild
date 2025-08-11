from intigers.lcg_v3 import * 
import unittest

def LCGV3_sample_k_(prg):
    start = 515 
    m = 311 
    a = -422 
    l0 = -299 
    l1 = 478 
    sc_size = 100 
    preproc_gd = False 
    super_range = [-300,500]
    return LCGV3(start,m,a,l0,l1,\
        sc_size,preproc_gd,prg,\
        super_range) 

def LCGV3_sample_k(): 
    prg = prg__LCG(17,59,-43,900) 
    return LCGV3_sample_k_(prg) 

def LCGV3_sample_q(): 
    prg = prg__LCG(76,599,-433,900) 
    return LCGV3_sample_k_(prg) 

### lone file test 
"""
python -m tests.test_lcg_v3 
"""
###
class LCGV3Methods(unittest.TestCase):

    def test__modrange_for_congruence__case1(self):

        modrange = [3,37]
        qx = modrange_for_congruence(15,25,modrange)
        assert modulo_in_range(25,qx) == 15 

        qx2 = modrange_for_congruence(15,5,modrange)
        assert modulo_in_range(5,qx2) == 15 
        
        qx3 = modrange_for_congruence(15,54,modrange)
        assert modulo_in_range(54,qx3) == 15 

        qx4 = modrange_for_congruence(15,-54,modrange)
        assert modulo_in_range(-54,qx4) == 15 

        modrange = [-30,37]
        qx5 = modrange_for_congruence(-15,-54,modrange)
        assert modulo_in_range(-54,qx5) == -15 

        qx6 = modrange_for_congruence(15,-549,modrange)
        assert modulo_in_range(-549,qx6) == 15 

        modrange = [-18321,143451]
        qx7 = modrange_for_congruence(-1500,-564321,modrange)
        assert modulo_in_range(-564321,qx7) == -1500 

        qx8 = modrange_for_congruence(1500,-564321,modrange)
        assert modulo_in_range(-564321,qx8) == 1500 

        qx9 = modrange_for_congruence(-1500,564321,modrange)
        assert modulo_in_range(564321,qx9) == -1500 

        qx10 = modrange_for_congruence(1500,111321,modrange)
        assert modulo_in_range(111321,qx10) == 1500 

        qx11 = modrange_for_congruence(-1500,111321,modrange)
        assert modulo_in_range(111321,qx11) == -1500 
        return  

    def test__unimodular_number__modulo_range_adjustment__case1(self):
        # case 1 
        wv = 38 
        av = 29 
        cv = 32 
        modrange = [3,38]
        mrx4 = unimodular_number__modulo_range_adjustment(wv,av,cv,modrange)
        assert modulo_in_range(29,mrx4) == wv 

        # case 2 
        wv = 5 
        av = 29 
        cv = 32 
        modrange = [3,38]
        mrx5 = unimodular_number__modulo_range_adjustment(wv,av,cv,modrange)
        assert modulo_in_range(29,mrx5) == wv 

        # case 3 
        wv = 13 
        cv = -12 
        av = 85 
        mrx2_ = unimodular_number__modulo_range_adjustment(wv,cv,av,[50,97])
        assert modulo_in_range(cv,mrx2_) == 13 

        # case 4 
        wv = 100 
        cv = -20 
        av = 77 
        mrx3_ = unimodular_number__modulo_range_adjustment(wv,cv,av,[50,97])
        assert modulo_in_range(-20,mrx3_) == 100 

        # case 5 
        wv = 50 
        cv = -20 
        av = 77 
        mrx4_ = unimodular_number__modulo_range_adjustment(wv,cv,av,[50,97])
        assert modulo_in_range(-20,mrx4_) == 50 

        return

    def test__multimodular_number__modulo_range_adjustment__case1(self):
        
        # case 1
        modrange = [3,38]
        pv = 7 
        av = 5000 # 33 
        sign = -1 
        mrx6 = multimodular_number__modulo_range_adjustment(pv,av,sign,modrange)
        assert to_trinary_relation(modulo_in_range(av,mrx6),pv) == sign 

        # case 2 
        pv = 37
        sign = 1
        mrx7 = multimodular_number__modulo_range_adjustment(pv,av,sign,modrange)
        assert to_trinary_relation(modulo_in_range(av,mrx7),pv) == sign

        # case 3 
        pv = 15 
        av = 5010 
        sign = 1 
        mrx8 = multimodular_number__modulo_range_adjustment(pv,av,sign,modrange)
        assert to_trinary_relation(modulo_in_range(av,mrx8),pv) == sign 
 
    def test__LCGV3__adjust_modulo_range__case1(self): 
        lcg1 = prg__LCG(13,-20,14,1000)
        lg = LCGV3(7,10,5,3,38,10)

        # case 1 
        rv = 8 
        nv = 36 
        sign = -1 
        mr23 = lg.adjust_modulo_range(rv,nv,sign,lcg1)
        avx = lg.m * rv + lg.a
        assert to_trinary_relation(modulo_in_range(avx,mr23),rv) == sign 

        # case 2 
        rv = 16 
        sign = -1 
        mr23 = lg.adjust_modulo_range(rv,42,sign,lcg1)
        avx = lg.m * rv + lg.a
        assert to_trinary_relation(modulo_in_range(avx,mr23),rv) == sign 

        # case 3 
        lg2 = LCGV3(7,-2,6,3,38,10)
        sign = -1 
        mr24 = lg2.adjust_modulo_range(4,36,sign,lcg1)
        avx = lg2.m * rv + lg2.a
        assert to_trinary_relation(modulo_in_range(avx,mr24),rv) == sign 

        # case 4 
        lg3 = LCGV3(7,-1,10,3,38,10)
        rv = 6 
        sign = 1 
        mr25 = lg3.adjust_modulo_range(rv,4,sign,lcg1)
        avx = lg3.m * rv + lg3.a
        assert to_trinary_relation(modulo_in_range(avx,mr25),rv) == sign 

        # case 6 
        sign = 0 
        mr26 = lg3.adjust_modulo_range(rv,4,sign,lcg1)
        avx = lg3.m * rv + lg3.a
        assert to_trinary_relation(modulo_in_range(avx,mr26),rv) == sign 

        # case 7 
        sign = 0 
        rv = 36
        avx = lg3.m * rv + lg3.a
        mr27 = lg3.adjust_modulo_range(rv,12,sign,lcg1)
        assert to_trinary_relation(modulo_in_range(avx,mr27),rv) == sign 

    def test__LCGV3__next__case1(self): 
        lcg1 = prg__LCG(13,-20,14,1000)
        lg2 = LCGV3(7,-2,6,3,38,10)
        lg2.set_prg(lcg1)
        lg2.autoset_tv(1,6,ext_prg=lcg1) 

        lx = [] 
        for _ in range(7): 
            lx_ = next(lg2) 
            lx.append(lx_)

        lx = np.array(lx) 
        sv = stdop_vec(lx,to_trinary_relation,cast_type=np.int32)
        assert lg2.tv == sv 

    def test__LCGV3__next__case2(self): 

        g3 = LCGV3_sample_k() 
        gen_type = 2 
        l = 7 
        is_delta2_mutable = True 
        g3.static_autoset(gen_type,l,is_delta2_mutable)

        qs = [] 
        cx = 0 

        stat = True 
        while stat: 
            q = next(g3) 
            qs.append(q) 
            cx += bool(g3.stat__new_trinary) 
            stat = not g3.stat__new_trinary 

        assert len(qs) == 58, "got {}".format(len(qs))

    """
    ensures that LCGV3 changes its trinary vector at the 
    appropriate rate. 
    """
    def test__LCGV3__next_batch__auto_td__case1(self): 
        g3 = LCGV3_sample_q() 
        gen_type = 2 
        l = 7 
        is_delta2_mutable = True 
        g3.static_autoset(gen_type,l,is_delta2_mutable)

        d0,d1 = g3.delta_one,g3.delta_two

        bx0 = g3.next_batch__auto_td()
        bx1 = g3.next_batch__auto_td()
        bx2 = g3.next_batch__auto_td()
        bx3 = g3.next_batch__auto_td()
        bx4 = g3.next_batch__auto_td()

        d0_,d1_ = g3.delta_one,g3.delta_two

        assert (d0,d1) == (7,4) and (d0_,d1_) == (4,7)

if __name__ == '__main__':
    unittest.main()