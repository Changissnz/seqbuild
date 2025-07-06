from desi.not_equals import * 
from morebs2.numerical_generator import prg__LCG,prg__n_ary_alternator

import unittest

### lone file test 
"""
python3 -m tests.test_not_equals
"""
###
class NotEqualsMethods(unittest.TestCase):

    def test__not_equals__pairvec__case1(self):

        lcg_px = prg__LCG(12,3,13,250)

        v1 = np.array([3.0,13.4,26.5])
        v2 = np.array([3.0,113.4,26.5])

        nep = not_equals__pairvec(v1,v2,lcg_px,indices=None)
        assert not np.any(nep[0] == nep[1])

        lcg_px2 = prg__LCG(321,43,117,500)

        cx = set() 
        for _ in range(50): 
            nep2 = not_equals__pairvec(v1,v2,lcg_px2) 
            assert not np.any(nep2[0] == nep2[1])
            cx |= {vector_to_string(nep2[1],float)}

        assert len(cx) == 20 
        assert cx == \
        {'4.009375,114.75437500000001,27.5828125',\
            '4.0108303249097474,114.80938628158846,27.595667870036102', \
            '4.0062893081761,114.63773584905661,27.555555555555557', \
            '4.007125890736342,114.66935866983374,27.56294536817102', \
            '4.024793388429752,115.33719008264464,27.71900826446281', \
            '4.038961038961039,115.87272727272727,27.844155844155843', \
            '4.107142857142858,118.45,28.446428571428573', \
            '4.142857142857142,119.80000000000001,28.761904761904763', \
            '4.007957559681698,114.70079575596817,27.57029177718833', \
            '4.025,115.345,27.720833333333335', \
            '4.0131578947368425,114.89736842105263,27.61622807017544', \
            '4.0234375,115.2859375,27.70703125', \
            '4.013574660633484,114.9131221719457,27.619909502262445', \
            '4.007142857142857,114.67,27.563095238095237', \
            '4.013636363636364,114.91545454545455,27.620454545454546', \
            '4.016949152542373,115.0406779661017,27.649717514124294', \
            '4.007009345794392,114.66495327102804,27.561915887850468', \
            '4.009146341463415,114.74573170731708,27.58079268292683', \
            '4.009345794392523,114.75327102803739,27.582554517133957', \
            '4.15,120.07000000000001,28.825'}
        

        lcg_px3 = prg__n_ary_alternator(3,5000,5)

        cx2 = set() 
        for _ in range(50): 
            nep3 = not_equals__pairvec(v1,v2,lcg_px3) 
            assert not np.any(nep3[0] == nep3[1])
            cx2 |= {vector_to_string(nep3[1],float)}

        assert len(cx2) == 50 

        return 
    
    def test__not_equals__matrix_whole__case1(self):
        M = np.ones((5,3),dtype=float)
        M = M * 30.2

        prg = prg__LCG(-2,131,40,4001)
        submat_type = "L+L"

        qm = not_equals__matrix_whole(deepcopy(M),prg,submat_type)

        submat_type = "L+U"
        qm2 = not_equals__matrix_whole(deepcopy(M),prg,submat_type)

        submat_type = "R+L"
        qm3 = not_equals__matrix_whole(deepcopy(M),prg,submat_type)

        submat_type = "R+U"
        qm4 = not_equals__matrix_whole(deepcopy(M),prg,submat_type)

        sol = np.array([[31.1 , 30.2 , 30.2 ],\
            [32.  , 30.2 , 30.2 ],\
            [33.05, 30.5 , 30.2 ],\
            [35.  , 32.3 , 32.45],\
            [38.45, 34.85, 34.7 ]]) 

        assert np.all(np.round(np.abs(sol - qm2),5) == 0.0)


        assert not np.any(qm == qm2) 
        assert not np.any(qm == qm3) 
        assert not np.any(qm == qm4) 
        assert not np.any(qm2 == qm3) 
        assert not np.any(qm2 == qm4) 

        l = np.where(qm3 == qm4)
        assert len(l[0]) == 3 


        M2 = np.ones((1,11),dtype=float)
        M2 *= 40.7

        for x in ["L+L","L+U","R+L","R+U"]:
            qx = not_equals__matrix_whole(deepcopy(M2),prg,x)
            assert len(np.unique(qx)) == len(qx.flatten()) 

    def test__io_closest_multiples_MULO(self): 
        x = 300.0 
        y = 501.0 
        # case 1 
        q = io_closest_multiples_MULO(x,y,unit_epsilon=10**-2)
        assert q == (1.66,1.68)

        # case 2 
        q2 = io_closest_multiples_MULO(x,y,unit_epsilon=10**-1)
        assert y >= q2[0] * x and y <= q2[1] * x 

        # case 3 
        q3 = io_closest_multiples_MULO(x,y,unit_epsilon=1.0)
        assert q3 == (1,2)

        # case 4 
        q4 = io_closest_multiples_MULO(x,-y,unit_epsilon=1.0)
        assert q4 == (-2,-1)

        # case 5
        x2 = 6.0 
        y2 = np.array([4,14,50,66.0])
        q5 = io_closest_multiples_MULO(x2,y2,unit_epsilon=1.0)

        b = np.array([x2 * q5[0],x2 * q5[1]]).T 
        assert point_in_bounds(b,y2)
        qb = b[:,1] - b[:,0] 
        qd = safe_div(qb,1.0) 
        assert np.all(np.array(qd,dtype=int) == qd)

        # case 6
        q6 = io_closest_multiples_MULO(x2,y2,unit_epsilon=2.0)
        b = np.array([x2 * q6[0],x2 * q6[1]]).T 
        assert np.all(b[:,1] - b[:,0] == 6.0 * 2.0) 
        qb = b[:,1] - b[:,0] 
        qd = safe_div(qb,2.0) 
        assert np.all(np.array(qd,dtype=int) == qd)

        # case 7 
        q7 = io_closest_multiples_MULO(y2,x2,unit_epsilon=1.0)
        b = np.array([y2 * q7[0],y2 * q7[1]]).T 
        xx = np.ones((4,)) * x2 
        assert point_in_bounds(b,xx)

        # case 8 
        y2[0] *= -1 
        y2[2] *= -1 
        q8 = io_closest_multiples_MULO(y2,x2,unit_epsilon=1.0)
        b = np.array([y2 * q8[0],y2 * q8[1]]).T 
        b = np.sort(b,axis = 1)
        assert point_in_bounds(b,xx) 

        # case 9 
        b2 = io_closest_range_MULO(y2,x2,unit_epsilon=1.0) 
        assert equal_iterables(b,b2)

if __name__ == '__main__':
    unittest.main()