from mini_dm.vs_fit import * 
from morebs2.numerical_generator import prg__LCG, prg__constant,prg__n_ary_alternator
import unittest

### lone file test 
"""
python -m tests.test_vs_fit
"""
###
class VSFitMethods(unittest.TestCase):


    def test__ratio__type_asymmetricANDsymmetric(self): 

        q0,q1 = 3,4.0
        qx = ratio__type_asymmetric(q0,q1,"min 1.0")
        assert qx == 4 / 3
        qx2 = ratio__type_asymmetric(q0,q1,"max 1.0")
        assert qx2 == 3 / 4

        qx3 = ratio__type_symmetric(q0,q1,ref=0)
        qx4 = ratio__type_symmetric(q0,q1,ref=1)
        assert qx3 + qx4 == 1.0 

        q0,q1 = 0.,2. 
        qx5 = ratio__type_symmetric(q0,q1,ref=1)
        qx6 = ratio__type_symmetric(q0,q1,ref=0)
        assert qx5 == 1.0 and qx6 == 0.0 

        qx7 = ratio__type_asymmetric(q0,q1,"min 1.0")
        qx8 = ratio__type_asymmetric(q0,q1,"max 1.0")
        assert qx7 == qx8 and qx7 == 0.0 

    def test__ratio_vector__case1(self): 

        q0 = np.array([4.0,2.1,7.7])
        q1 = np.array([5.0,6.3,10.90])
        rtype = "a" 
        parameter = "min 1.0" 
        parameter2 = -1 
        rv1 = ratio_vector(q0,q1,rtype,parameter,parameter2)

        parameter2 = 0 
        rv2 = ratio_vector(q0,q1,rtype,parameter,parameter2)

        parameter2 = 1
        rv3 = ratio_vector(q0,q1,rtype,parameter,parameter2)

        rv_sol = np.array([1.25      , 3.        , 1.41558442])
        assert np.all(abs(rv1 - rv_sol) <= 10**-5)

        assert set(np.unique(rv2)) == {1.0} 
        assert np.all(rv1 == rv3) 


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