from mini_dm.vs_fit_mem import * 
from morebs2.matrix_methods import equal_iterables
import unittest


### lone file test 
"""
python -m tests.test_vs_fit_mem
"""
###
class VSFitMemMethods(unittest.TestCase):

    def test__MAHypMach__io_to_AffineDelta__case1(self):
        # case 1 
        a2 = np.array([500,1000,1000,1120,2500])
        mx = np.array([100,4200,2000,200,-700])
        ad6 = AffineDelta(mx,a2,0) 

        x = 45 
        y = ad6.fit(x) 
        d = ad6.abssum_diff(x) 
        cv = ad6.cvec(x) 

        ma_dim = [5,5] 
        adx = MAHypMach.io_to_AffineDelta(x,y,d,cv,ma_dim,ma_order=0)

        assert np.all(adx.fit(x) == y) 

        # case 2 
        a2 = np.array([500,1000,1000,1120,-2500])
        mx = np.array([100,4200,2000,200,700])
        ad7 = AffineDelta(mx,a2,0) 
        cv = ad7.cvec(x) 
        y = ad7.fit(x) 
        d = ad7.abssum_diff(x) 
        ma_order = 0 

        adx3 = MAHypMach.io_to_AffineDelta(x,y,d,cv,ma_dim,ma_order)
        assert np.all((adx3.fit(45) == adx.fit(45)) == \
            np.array([ True,  True,  True,  True, False]))

        # case 3 
        a2 = np.array([-50,1000,1000,1120,-2500])
        mx = np.array([-19,420,-2000,200,700])
        ad8 = AffineDelta(mx,a2,0) 
        cv = ad8.cvec(x) 
        y = ad8.fit(x) 
        d = ad8.abssum_diff(x) 
        ma_order = 0 

        adx4 = MAHypMach.io_to_AffineDelta(x,y,d,cv,ma_dim,ma_order)
        assert equal_iterables(adx4.m,ad8.m) 
        assert equal_iterables(adx4.a,ad8.a) 

        # case 4 
        ad9 = AffineDelta(mx,a2,1)
        cv = ad9.cvec(x) 
        y = ad9.fit(x) 
        d = ad9.abssum_diff(x) 
        ma_order = 1
        adx5 = MAHypMach.io_to_AffineDelta(x,y,d,cv,ma_dim,ma_order)
        assert np.all(adx5.m == ad9.m)
        assert np.all(adx5.a == ad9.a)

        # case 5 
        a2 = 50 
        mx = -19 
        ad10 = AffineDelta(mx,a2,0) 
        cv = ad10.cvec(x) 
        y = ad10.fit(x) 
        d = ad10.abssum_diff(x) 
        ma_order = 0 

        adx6 = MAHypMach.io_to_AffineDelta(x,y,d,cv,[0,0],ma_order)
        assert adx6.m == mx 
        assert adx6.a == a2 

        # case 6 
        a2 = 0 
        mx = -19 
        ad11 = AffineDelta(mx,a2,0) 
        cv = ad11.cvec(x) 
        y = ad11.fit(x) 
        d = ad11.abssum_diff(x) 
        ma_order = 0 

        adx7 = MAHypMach.io_to_AffineDelta(x,y,d,cv,[0,0],ma_order)
        assert adx7.a == a2 
        assert adx7.m == mx 

        # case 7 
        a2 = 19 
        mx = 0  
        ad12 = AffineDelta(mx,a2,0) 
        cv = ad12.cvec(x) 
        y = ad12.fit(x) 
        d = ad12.abssum_diff(x) 
        ma_order = 0 

        adx8 = MAHypMach.io_to_AffineDelta(x,y,d,cv,[0,0],ma_order)
        assert adx8.a == a2 
        assert adx8.m == mx 
        return 

if __name__ == '__main__':
    unittest.main()

