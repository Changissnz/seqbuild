from intigers.lcg_v2 import * 
import unittest

### lone file test 
"""
python -m tests.test_lcg_v2 
"""
###
class LCGV2Methods(unittest.TestCase):

    def test__LCGV2__basefunctions__case1(self):

        # case 1
        l2 = LCGV2(2,2,0,0,10,0) 

        d1 = l2.init_dir()
        d2 = l2.gen_dir() 
        assert d1 == 1 
        assert d2 == 1 

        cv = l2.coverage()
        assert cv == 0.4 

        cl = l2.cycle_length()
        assert cl == 4 

        # case 2 
        l2 = LCGV2(0,2,1,0,10,0) 

        d1 = l2.init_dir()
        d2 = l2.gen_dir() 
        assert d1 == 1 
        assert d2 == 1 

        cv = l2.coverage()
        assert cv == 0.45 

        cl = l2.cycle_length()
        assert cl == 5 


if __name__ == '__main__':
    unittest.main()