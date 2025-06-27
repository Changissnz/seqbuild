from mini_dm.vs_fit import * 
from morebs2.numerical_generator import prg__LCG, prg__constant,prg__n_ary_alternator
import unittest

### lone file test 
"""
python -m tests.test_vs_fit
"""
###
class VSFitMethods(unittest.TestCase):

    def test__AffineDelta__next__case1(self):

        m_type = "vec"
        a_type = "vec" 
        prg = prg__LCG(71,688,31,900) 
        r_out1 = prg__constant((50.,5000.))
        r_out2 = prg__constant((50.,5000.))
        ma_order = 0

        ad = AffineDelta.one_instance(m_type,a_type,prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)

        sol = "m: [329. 333. 385. 161.]\na: [849. 793.  65. 501.]\no: 0\n"
        assert str(ad) == sol 

        ad2 = AffineDelta.one_instance_(prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)

        prg = prg__constant(1) 
        prg = prg__n_ary_alternator(3,500,5)
        ad3 = AffineDelta.one_instance_(prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)
        ad4 = AffineDelta.one_instance_(prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)

        prg = prg__constant(2) 
        ad5 = AffineDelta.one_instance_(prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)

        assert ad.type() == (1, 1)
        assert ad2.type() == (0, 0)
        assert ad3.type() == (1, 0)
        assert ad4.type() == (0, 1)
        assert ad5.type() == (1, 1)


if __name__ == '__main__':
    unittest.main()