from intigers.lcg_v3 import * 
import unittest

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


if __name__ == '__main__':
    unittest.main()