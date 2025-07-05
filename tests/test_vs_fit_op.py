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
        #assert ad_diff.m != 0.0 
        assert not np.any(ad_diff.m == 0.)
        assert not np.any(ad_diff.a == 0)

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

    def test__IOFit__init_MAHypMach__case1(self):
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
        dx = np.array([vst.ad.abssum_diff(x_) for x_ in X]) 
        cv = [ad.cvec(x_) for x_ in X]
        ma_dim = [0,4]
        ma_order = 0 
        iof.load_mahyp_auxvar(dx,cv,ma_dim,ma_order)

        iof.init_MAHypMach()

        assert iof.mahm.mhm.indices == [0,1,2,3,4]
        ix = iof.mahm.mhm.info[0]
        for x in iof.mahm.mhm.info: 
            assert x == ix 
            print(x)

    def test__IOFit__io_pointpair_to_AffineDelta__case1(self):

        m_type = "vec"
        a_type = "vec" 
        prg = prg__LCG(71,688,31,900) 
        r_out1 = prg__constant((50.,5000.))
        r_out2 = prg__constant((50.,5000.))
        ma_order = 0

        ad = AffineDelta.one_instance(m_type,a_type,prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)

        x1,x2 = 56,400 
        y1 = ad.fit(x1)
        y2 = ad.fit(x2)

        ad2 = IOFit.io_pointpair_to_AffineDelta(x1,x2,y1,y2,ma_order=0)

        ad.ma_order = 1
        y1 = ad.fit(x1) 
        y2 = ad.fit(x2) 
        ad3 = IOFit.io_pointpair_to_AffineDelta(x1,x2,y1,y2,ma_order=1)

        ad2.ma_order = 1 
        assert ad == ad2 
        assert ad2 == ad3 
        return

if __name__ == '__main__':
    unittest.main()