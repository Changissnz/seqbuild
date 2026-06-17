from intigers.modulo_ops import * 
import unittest

### lone file test 
"""
py -m tests.test_modulo_ops
""" 
###
class ModuloOpsMethods(unittest.TestCase):
 
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
        assert modulo_in_range(cv,mrx2_) == 13, "got {}".format(modulo_in_range(cv,mrx2_))

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
 
 
if __name__ == '__main__':
    unittest.main()