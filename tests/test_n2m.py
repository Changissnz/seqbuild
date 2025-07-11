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

        qsol = np.array([-1,  0,  0, 0, 0])
        assert np.all(q == qsol)

        qsol2 = np.array([1,  0,  0, 0, 0])
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
        assert np.all(h_ == 0)
        assert np.all(h2_ == 0)

        # case 5
        ac.add(x0,x1,e0,e1) 
        g_ = ac.induce_derivative_v2(x1_,x0_)
        g2_ = ac.induce_derivative_v2(x2_,x0_)
        assert np.all(g_ == 0)
        assert np.all(g2_ == 0)

        ac.add(x0,x1,e0,e1) 
        ac.add(x0,x1,e0,e1) 
        ac.add(x0,x1,e0,e1) 
        ac.add(x0,x1,e0,e1) 

        # case 6
        j_ = ac.induce_derivative_v2(x1_,x0_)
        j2_ = ac.induce_derivative_v2(x2_,x0_)

        jsol1 = np.array([1, 0, 0, 0, 0])
        jsol2 = np.array([-1,  0,  0,0,0])

        assert np.all(j_ == jsol1)
        assert np.all(j2_ == jsol2)

    def test__N2MAutocorrelator__induce_derivative_case2(self):

        nm = (8,1) 
        ac = N2MAutocorrelator(nm) 

        x0 = np.array([5,1,41,0,0,1,1,1]) 
        x1 = np.array([10,2,45,-7,3,40,421,-34])  

        e0 = np.array([1]) 
        e1 = np.array([2])  
        ac.add(x0,x1,e0,e1)

        e0_,e1_= np.array([2]),np.array([1]) 
        for i in range(10): 
            ac.add(x0,x1,e0_,e1_)

        q = ac.induce_derivative(x0,x1)

        x2 = np.array([1,2,0,0,0,0,0,0])
        x3 = np.array([2,3,0,0,0,0,0,0])
        q2_ = ac.induce_derivative(x2,x3)  
        assert np.all(q2_ == np.array([0])) 

        x4 = np.array([1,2,1,1,0,0,0,0])
        x5 = np.array([2,3,0,0,0,0,0,0])
        q2__ = ac.induce_derivative(x4,x5)  
        assert np.all(q2__ == np.array([0])) 

        x6 = np.array([1,2,1,1,1,1,0,0])
        x7 = np.array([2,3,0,0,2,2,0,0])
        q2 = ac.induce_derivative(x6,x7)  

        for i in range(50): 
            ac.add(x0,x1,e0_,e1_)

        q3 = ac.induce_derivative(x6,x7)  
        q4 = ac.induce_derivative_v2(x6,x7)  
        assert np.all(q3 == 0)
        assert np.all(q4 == 0)

        e0_ = np.array([-1])
        e1_ = np.array([4])
        ac.add(x6,x7,e0_,e1_) 

        x8 = x6 + 2 
        x9 = x7 + 2
        q5 = ac.induce_derivative(x8,x9)   
        q6 = ac.induce_derivative_v2(x8,x9)   
        sol56 = np.array([1])
        assert np.all(q5 == sol56)
        assert np.all(q6 == sol56) 

        x10 = x6 + 1 
        x11 = x7 + 2
        q7 = ac.induce_derivative(x10,x11)   
        q8 = ac.induce_derivative_v2(x10,x11)   
        sol78 = np.array([0])
        assert np.all(q7 == sol78)
        assert np.all(q8 == sol78) 



if __name__ == '__main__':
    unittest.main()