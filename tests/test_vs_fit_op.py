from mini_dm.vs_fit_op import * 
from morebs2.numerical_generator import prg__LCG, prg__constant,prg__n_ary_alternator
from morebs2.matrix_methods import equal_iterables
import unittest


### lone file test 
"""
python -m tests.test_vs_fit_op
"""
###
class VSTransformMethods(unittest.TestCase):

    def test__VSTransform__diff_ad__case1(self):

        m4,a4 = np.array([440,1400,250,500]),np.array([20,30,40,50])
        ad4 = AffineDelta(m4,a4,0)

        m5,a5 = 0.05,0.05 
        ad5 = AffineDelta(m4+m5,a4+a5,0) 

        vst = VSTransform(ad4)

        q2 = vst.diff_ad(24,ad5,op_type="all",\
                dfunc=lambda x,x2: np.sum(np.abs(x - x2)))
        q2_ = vst.diff_ad(24,ad5,op_type="one",\
                dfunc=lambda x,x2: np.sum(np.abs(x - x2)))
        assert round(q2_ / q2) == 29,"got {}".format(q2_/q2)

        q3 = vst.diff_ad(45,ad5,op_type="all",\
                dfunc=lambda x,x2: np.sum(np.abs(x - x2)))
        q4 = vst.diff_ad(45,ad5,op_type="one",\
                dfunc=lambda x,x2: np.sum(np.abs(x - x2)))
        assert round(q4 / q3) == 16, "got {}".format(q4/q3)
        return
    
    def test__VSTransform__cmp_ad__case1(self):

        m4,a4 = np.array([440,1400,250,500]),np.array([20,30,40,50])
        ad4 = AffineDelta(m4,a4,0)

        m5,a5 = 0.05,0.05 
        ad5 = AffineDelta(m4+m5,a4+a5,0) 

        vst = VSTransform(ad4)
        vst2 = VSTransform(ad5)

        q = vst.cmp_ad(ad5)
        q2 = vst2.cmp_ad(ad4)

        assert q == q2 

        assert abs(round(q[0] - 0.4,5)) < 10 ** -5 
        assert abs(round(q[1] - 0.4,5)) < 10 ** -5 


    def test__IOFit__diff__case1(self):
        m,a = 34,np.array([4,-14.0,29,79.0])
        ad = AffineDelta(m,a,0)

        # case 1
        vst = VSTransform(ad) 
        hdf = vst.to_hypdiff_func() 
        maf = vst.to_ma_diff_func()

        X = np.array([45,52,61,15,150])
        Y = np.array([ad.fit(x) for x in X])

        unknownf = vst.ad.fit 

        p3 = {"s":0,"a":"max 1.0"}
        md = ad.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        md_hyp = MADHyp.naive_hyp_for_MADescriptor(md)

        ad2 = md_hyp.solve_into_AffineDelta([0,4],ma_order=0)
        vst2 = VSTransform(ad2) 
        maf = vst2.to_ma_diff_func() 


        iof = IOFit(X,Y,unknownf,hdf,maf)
        iof.load_hyp(md_hyp,[4,4])

        dx = set() 
        for x in X: 
            q = iof.hyp_diff(x) 
            dx |= {str(q)}

        assert len(dx) == len(X)  

        ad_diff = iof.ma_diff(ad)
        assert ad_diff.m != 0.0 
        assert not np.any(ad_diff == 0.)

        # case 2 
        X20 = [4,5,6,7]
        X21 = [-3,13,10,9] 
        X2 = np.array([X20,X21])
        Y2 = np.array([ad.fit(x) for x in X2])

        iof2 = IOFit(X2,Y2,unknownf,hdf,maf)
        iof2.load_hyp(md_hyp,[4,4])

        S2 = set() 
        for x in X2: 
            q = iof.hyp_diff(x) 
            assert len(q) == 4 
            S2 |= {str(q)}
        assert len(S2) == len(X2)

if __name__ == '__main__':
    unittest.main()