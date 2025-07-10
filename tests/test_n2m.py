from mini_dm.n2m import * 
import unittest

### lone file test 
"""
python -m tests.test_n2m
"""
###
class N2MAutocorrelatorMethods(unittest.TestCase):

    def test__N2MAutocorrelator__induce_derivative_case1(self):

        nm = (5,5)
        ac = N2MAutocorrelator(nm)


        x0 = np.array([5,1,41,0,0]) 
        x1 = np.array([10,2,45,0,0]) 

        e0 = np.array([1,1,1,5,5]) 
        e1 = np.array([2,1,1,4,2]) 

        ac.add(x0,x1,e0,e1)
        ac.add(x0,x1,e0,e1)
        ac.add(x0,x1,e0,e1)
        ac.add(x0,x1,e0,e1)
        ac.add(x0,x1,e0,e1)

        x0_ = np.array([1,0,0,0,0])
        x1_ = np.array([2,0,0,0,0])
        x2_ = np.array([0,0,0,0,0])

        # case 1 
        q = ac.induce_derivative(x1_,x0_)
        q2 = ac.induce_derivative(x2_,x0_)

        qsol = np.array([-1,  0,  0, -1, -1])
        assert np.all(q == qsol)

        qsol2 = np.array([1,  0,  0, -1, -1])
        assert np.all(q2 == qsol2)

        # case 2 
        q_ = ac.induce_derivative_v2(x1_,x0_)
        q2_ = ac.induce_derivative_v2(x2_,x0_)
        assert np.all(q_ == qsol)
        assert np.all(q2_ == qsol2)

        x0 = np.array([-10,0,0,0,0])
        x1 = np.array([1,0,0,0,0]) 
        e0 = np.array([2,0,0,0,0])
        e1 = np.array([1,0,0,0,0])

        ac.add(x0,x1,e0,e1) 

        # case 3 
        h = ac.induce_derivative(x1_,x0_) 
        h2 = ac.induce_derivative(x2_,x0_) 
        hsol = np.array([1, 0, 0, 0, 0])
        hsol2 = np.array([-1,  0,  0,  0,  0])
        assert np.all(h == hsol) 
        assert np.all(h2 == hsol2) 

        # case 4 
        h_ = ac.induce_derivative_v2(x1_,x0_)
        h2_ = ac.induce_derivative_v2(x2_,x0_)
        hsol_ = np.array([0,  0,  0, -1, -1])
        hsol2_ = np.array([0,  0,  0, -1, -1])
        assert np.all(h_ == hsol_)
        assert np.all(h2_ == hsol2_)

        # case 5
        ac.add(x0,x1,e0,e1) 
        g_ = ac.induce_derivative_v2(x1_,x0_)
        g2_ = ac.induce_derivative_v2(x2_,x0_)
        assert np.all(g_ == hsol2_)
        assert np.all(g2_ == hsol2_)

        ac.add(x0,x1,e0,e1) 
        ac.add(x0,x1,e0,e1) 
        ac.add(x0,x1,e0,e1) 
        ac.add(x0,x1,e0,e1) 

        # case 6 
        j_ = ac.induce_derivative_v2(x1_,x0_)
        j2_ = ac.induce_derivative_v2(x2_,x0_)

        jsol1 = np.array([1, 0, 0, 0, 0])
        jsol2 = np.array([-1, 0, 0,0,0])

        assert np.all(j_ == jsol1)
        assert np.all(j2_ == jsol2)

if __name__ == '__main__':
    unittest.main()