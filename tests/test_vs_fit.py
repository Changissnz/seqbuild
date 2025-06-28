from mini_dm.vs_fit import * 
from morebs2.numerical_generator import prg__LCG, prg__constant,prg__n_ary_alternator
from morebs2.matrix_methods import equal_iterables
import unittest

def AffineDelta__sampleA():
    m_type = "vec"
    a_type = "vec" 
    prg = prg__LCG(71,688,31,900) 
    r_out1 = prg__constant((50.,5000.))
    r_out2 = prg__constant((50.,5000.))
    ma_order = 0
    ad = AffineDelta.one_instance(\
        m_type,a_type,prg,r_out1,r_out2,dim_range=None,\
        ma_order=ma_order)
    return ad 

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

    def test__ratio_vector__case2(self): 
        q0 = np.array([4.0,2.1,7.7])
        q1 = np.array([3.0,6.3,5.90])

        parameter = "min 1.0"
        parameter2 = -1 
        rv4 = ratio_vector(q0,q1,"a",parameter,parameter2)

        parameter2 = 0 
        rv5 = ratio_vector(q0,q1,"a",parameter,parameter2)

        parameter2 = 1 
        rv6 = ratio_vector(q0,q1,"a",parameter,parameter2)

        rx = np.array([rv5,rv6])

        parameter = "max 1.0"
        parameter2 = -1 
        rv7 = ratio_vector(q0,q1,"a",parameter,parameter2)

        parameter2 = 0 
        rv8 = ratio_vector(q0,q1,"a",parameter,parameter2)

        parameter2 = 1 
        rv9 = ratio_vector(q0,q1,"a",parameter,parameter2)

        rx2 = np.array([rv8,rv9])

        sol_rx = np.array([\
            [1.33333333, 1.        , 1.30508475],\
            [1.        , 3.        , 1.        ]])

        sol_rx2 = np.array([\
            [1.        , 0.33333333, 1.        ],\
            [0.75      , 1.        , 0.76623377]])

        assert equal_iterables(sol_rx,rx)
        assert equal_iterables(sol_rx2,rx2)

    # tests for rtype == "auto" 
    def test__ratio_vector__case3(self): 
        q0 = np.array([4.,5.,15.,-27.])
        q1 = np.array([14.,15.,9.,54.])
        rvp = {"a":"min 1.0","s":1}

        rvx1 = ratio_vector(q0,q1,"auto","min 1.0",0,rvp)

        rvp = {"a":"max 1.0","s":0}
        rvx2 = ratio_vector(q0,q1,"auto",None,0,rvp)

        rvp = {"a":"max 1.0","s":1}
        rvx3 = ratio_vector(q0,q1,"auto","min 1.0",0,rvp)

        rvp = {"a":"min 1.0","s":0}
        rvx4 = ratio_vector(q0,q1,"auto","min 1.0",0,rvp)

        rvp = {"a":"min 1.0","s":1}
        rvx5 = ratio_vector(q0,q1,"auto","min 1.0",0,rvp)

        sol2 = np.array([ 0.22222222,0.25,0.625,-0.5])
        sol3 = np.array([ 0.77777778,0.75,0.375,-0.5])
        sol4 = np.array([ 0.22222222,0.25,0.625,-2.])
        sol5 = np.array([ 0.77777778,0.75,0.375,-2.])

        assert equal_iterables(sol2,rvx2)
        assert equal_iterables(sol3,rvx3)
        assert equal_iterables(sol4,rvx4)
        assert equal_iterables(sol5,rvx5)

        # test w/ full output vector of element (float,type)
        rvx5 = ratio_vector(q0,q1,"auto","min 1.0",0,rvp,1) 

        sol6 = np.array([\
            ['0.7777777777777778', 's'],\
            ['0.75', 's'],\
            ['0.375', 's'],\
            ['-2.0', 'a']], dtype='<U32')

        assert np.all(sol6 == rvx5)
        return 

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

    def test__AffineDelta__cvec__case1(self):
        ad = AffineDelta__sampleA() 
        sol = "m: [329. 333. 385. 161.]\na: [849. 793.  65. 501.]\no: 0\n"
        assert str(ad) == sol 

        cvec = ad.cvec(12)
        cvec2 = ad.cvec(150)

        assert np.all(cvec[0] + cvec[1] - 1.0 <= 10** -5) 
        assert np.all(cvec2[0] + cvec2[1] - 1.0 <= 10** -5) 
        assert not np.any(cvec[0] == cvec2[0])
        assert not np.any(cvec[1] == cvec2[1])
        return 
    
    def test__AffineDelta__to_ma_descriptor__case1(self): 
        ad = AffineDelta__sampleA()
        p3 = {"s": 0,"a":"max 1.0"} 
        md = ad.to_ma_descriptor(p3)


        sol_md = "* RV: [0.27928693 0.29573712 0.85555556 0.24320242]\n" + \
                "* RVT: ['s' 's' 's' 's']\n* T: [-1 -1  1 -1]\n* S: [1 1 1 1]\n" +\
                "* D: [520. 460. 320. 340.]\n"
        assert sol_md == str(md )

        sol_ad = "m: [329. 333. 385. 161.]\n" + \
            "a: [849. 793.  65. 501.]\n" + \
            "o: 0\n"

        assert sol_ad == str(ad)

        md2 = ad.to_ma_descriptor(p3,d_operator=lambda x,x2: np.abs(x) + np.abs(x2)) 
        assert not np.any(md.d_vec == md2.d_vec)

    def test__MADescriptor__solve_case1(self): 
        # case 1
        x1,x2 = 3.0,7.0 
        ad = AffineDelta(x1,x2,0)

        p3 = {"a":"max 1.0","s":0}
        md = ad.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION) 

        assert md.solve([0,0]) == (x1,x2) 

        # case 2 
        ad2 = AffineDelta(x1,x2,1) 
        md2 = ad2.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)

        sol_md = "* RV: 0.3\n* RVT: s\n* T: -1\n* S: 1\n* D: 10.0\n"
        assert str(md) == sol_md

        sol_md2 = "* RV: 0.7\n* RVT: s\n* T: 1\n* S: 1\n* D: 10.0\n"
        assert str(md2) == sol_md2

        # case 3
        x1,x2 = 3.0,np.array([1.0,7.0,-5.0])
        ad3 = AffineDelta(x1,x2,0)

        md3 = ad3.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr = md3.solve([0,3])

        assert abs(qr[0] - x1) < 10 ** -5 
        assert equal_iterables(x2,qr[1])

        # case 4 
        x1,x2 = np.array([0.5,15.0,12.0]),np.array([1.0,7.0,-5.0])

        ad4 = AffineDelta(x1,x2,1) 

        md4 = ad4.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr2 = md4.solve([3,3]) 

        assert equal_iterables(qr2[0],x2),"want {},got {}".format(x2,qr2[0])
        assert equal_iterables(qr2[1],x1)

    def test__MADescriptor__solve_case2(self): 
        # case 1 
        x1,x2 = np.array([0.0,15.0,12.0]),np.array([1.0,0.0,-5.0])

        ad5 = AffineDelta(x1,x2,1) 
        p3 = {"a":"max 1.0","s":0}

        md5 = ad5.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr3 = md5.solve([3,3]) 

        assert equal_iterables(qr3[0],x2)
        assert equal_iterables(qr3[1],x1),"want {}, got {}".format(x1,qr3[1])

        # case 2 
        x1,x2 = np.array([-10.0,5.0,27.0]),15.0 

        ad6 = AffineDelta(x1,x2,1) 
        p3 = {"a":"max 1.0","s":0}

        md6 = ad6.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr4 = md6.solve([0,3]) 
        qr5 = md6.solve([3,0])

        assert qr4[0] == x2
        assert equal_iterables(qr4[1],x1) 

        assert np.all(np.round(qr5[0],5) == [15,15,15.0]) and qr5[1] == 27.0 

    def test__MADescriptor__solve_case3(self): 

        p3 = {"a":"max 1.0","s":0}
        x1,x2 = -15,0

        ad7 = AffineDelta(x1,x2,0) 
        md7 = ad7.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr6 = md7.solve([0,0])

        ad8 = AffineDelta(x1,x2,1)
        md8 = ad8.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr7 = md8.solve([0,0])

        assert qr6[0] == x1 and qr6[1] == x2 
        assert qr7[0] == x2 and qr7[1] == x1 

    def test__MADescriptor__solve_case4(self):
        
        p3 = {"a":"min 1.0","s":0}
        
        # case 1
        x1,x2 = 0.0,np.array([1.0,7.0,-5.0])
        ad3 = AffineDelta(x1,x2,0)

        md3 = ad3.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr2 = md3.solve([0,3]) 
        assert qr2[0] == 0.0
        assert equal_iterables(qr2[1],x2)

        # case 2 
        x1,x2 = 6.0,np.array([0.0,12.0,3.0])
        ad4 = AffineDelta(x1,x2,0)

        md4 = ad4.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr3 = md4.solve([0,3]) 
        assert qr3[0] == 6.0
        assert equal_iterables(qr3[1],x2)

        # case 3 

        x1,x2 = np.array([6.0,24,-2]),np.array([10.0,-12.0,3.0])
        ad4 = AffineDelta(x1,x2,0)

        md4 = ad4.to_ma_descriptor(p3,d_operator=DEFAULT_MA_DISTANCE_FUNCTION)
        qr3 = md4.solve([3,3]) 

        assert equal_iterables(qr3[0],x1)
        assert equal_iterables(qr3[1],x2)


if __name__ == '__main__':
    unittest.main()