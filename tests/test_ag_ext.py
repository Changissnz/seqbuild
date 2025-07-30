from mini_dm.ag_ext import * 
from .idf_sample_z import * 
import unittest

### lone file test 
"""
python -m tests.test_ag_ext 
"""
###
class APRNGGaugeV2Methods(unittest.TestCase):

    def test__APRNGGaugeV2__std_cat_entropy__case1(self):
        # case 1 
        l = [2,0,4,6,8,10]
        lq = iter(l)

        def f(): 
            return next(lq,None) 

        aprng = f
        frange = [0.0,11.0]
        pradius = 0.25
        ag = APRNGGaugeV2(aprng,frange,pradius)

        term_func = lambda x,x2: not type(x) == type(None)
        mc = ag.measure_cycle(max_size=1000,term_func=term_func)

        is1 = IntSeq(ag.cycle)
        qx = ag.std_cat_entropy(is1,seg_length=2,start_value=None,\
                count_type="absdiff")
        qx2 = ag.std_cat_entropy(is1,seg_length=2,start_value=None,\
                count_type="equals")

        assert ag.catvec == [1, 0, 2, 3, 4, 5]
        assert mc == [0.25, 0.42424]
        assert qx == 0.24 
        assert qx2 == 1.0 

        # case 2 
        l = [1,3]
        lq = iter(l)
        def g(): 
            return next(lq,None) 

        ag.reload_var("aprng",g)
        ag.reload_var("frange",[0,4]) 
        ag.reload_var("pradius",1.0)

        ag.cycle = None 
        mc2 = ag.measure_cycle(max_size=1000,term_func=term_func)

        is2 = IntSeq(ag.cycle)
        qx3 = ag.std_cat_entropy(is2,seg_length=None,start_value=None,\
                count_type="absdiff")
        qx4 = ag.std_cat_entropy(is2,seg_length=None,start_value=None,\
                count_type="equals")

        assert mc2 == [1.0,0.5]
        assert qx3 == 1.0
        assert qx4 == 1.0
        assert ag.catvec == [0, 1]

        # case 3 
        l = [2,4,8,16,8,4,2,4,8,16]
        lq = iter(l)
        def h(): 
            return next(lq,None) 

        ag.reload_var("aprng",h)
        ag.cycle = None 
        mc3 = ag.measure_cycle(max_size=1000,term_func=term_func)
        is3 = IntSeq(ag.cycle)
        qx3 = ag.std_cat_entropy(is3,seg_length=None,start_value=None,\
                count_type="absdiff")
        assert ag.catvec == [0, 1, 2, 3, 2, 1, 0, 1, 2, 3]
        return

    def test__APRNGGaugeV2__std_cat_entropy__case2(self):
        L = [i for i in range(2,100)]
        prg1 = prg__LCG(13,87,211,5000)
        prg3 = prg__LCG(31,83,311,1500)

        def f(): 
            return modulo_in_range(prg1() * prg3() + prg1(),[0,5000]) 
        prg4 = f
        idf,_ = IDecForest_sampleZ(L,False,prg4)
        idf.verbose = False 

        def g(): 
            return next(idf) 

        aprng = g
        frange = [0.0,11.0]
        pradius = 0.25 
        ag = APRNGGaugeV2(aprng,frange,pradius)

        term_func = lambda x,x2: not type(x) == type(None)
        mc = ag.measure_cycle(max_size=100,term_func=term_func)

        is1 = IntSeq(ag.cycle)
        qx = ag.std_cat_entropy(is1,seg_length=None,start_value=None,\
                count_type="absdiff")
        assert max(ag.catvec) == 3 

    def test__APRNGGaugeV2__match_two_intseq(self):
        i1 = IntSeq([2,4,6,10])
        i2 = IntSeq([3,6,20])
        mf = absdiff_match_func
        match = APRNGGaugeV2.match_two_intseq(i1,i2,match_func=mf)

        assert match[2] == [3] 
        assert match[4] == [3] 
        assert match[6] == [6] 
        assert match[10] == [6] 
        q1 = APRNGGaugeV2.measure_match(match)
        assert q1[0] == 0 
        assert q1[1] == [6,3]

        i3 = IntSeq([1,3,4,8,9,11])
        mf = absdiff_match_func
        match = APRNGGaugeV2.match_two_intseq(i1,i3,match_func=mf)
        assert match[2] == [1,3] 
        assert match[4] == [4] 
        assert match[6] == [4,8] 
        assert match[10] == [9,11] 

        q2 = APRNGGaugeV2.measure_match(match)
        assert q2[0] == 3 
        assert q2[1] == [4]

    def test__APRNGGaugeV2__measure_matrix_chunk__case1(self):
        # case 1 
        d1 = 3 
        d0 = 5 
        m = None 

        m_ = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
        prg = prg__iterable(m_)
        ag = APRNGGaugeV2(prg,frange=(0.,1.),pradius=0.5) 
        d = ag.measure_matrix_chunk(m,d0,d1,axes={0,1}) 
        dsol00 = ([np.float64(1.0), np.float64(0.66667)], \
            0.5, (np.int32(1), np.int32(2), np.float64(1.3333333333333333)))
        assert np.all(np.round(np.array(d[0][0][0]) - dsol00[0]) <= 10 ** -5)
        assert d[0][0][1] == dsol00[1] 
        assert np.all(np.round(np.array(d[0][0][2]) - dsol00[2]) <= 10 ** -5)

        # case 2 
        e1 = [np.float64(250296.0533854167), np.float64(281979.40234375), \
            np.float64(283329.5872395834), np.float64(283616.1731770834), \
            np.float64(284118.4674479166), np.float64(284310.3658854166), \
            np.float64(284608.9049479166), np.float64(286900.2955729167), \
            np.float64(287596.70703125), np.float64(290375.7799479167)]
        e2 = [np.float64(328849.97374206624), np.float64(328852.83832539956), \
            np.float64(328854.1404087329), np.float64(328857.78624206624), \
            np.float64(328861.6924920663), np.float64(328863.2549920663), \
            np.float64(328863.5154087329), np.float64(328863.77582539956), \
            np.float64(328864.55707539956), np.float64(328866.6404087329)]
        m_ = e1 + e2 
        prg2 = prg__iterable(m_)
        ag2 = APRNGGaugeV2(prg2,frange=(0.,1.),pradius=0.5) 
        d2 = ag2.measure_matrix_chunk(None,2,len(e1),axes={0,1}) 
        d2sol10 = ([np.float64(1.0), np.float64(1.0)], \
            1.0, (np.int32(78553), np.int32(78553), np.float64(78553.0)))
        assert np.all(np.round(np.array(d2[1][0][0]) - d2sol10[0]) <= 10 ** -5)
        assert d2[1][0][1] == d2sol10[1] 
        assert np.all(np.round(np.array(d2[1][0][2]) - d2sol10[2]) <= 10 ** -5)

class AGExtOtherMethods(unittest.TestCase):

    def test__ranged_delta_decomposition__case1(self): 

        S = np.array([13,14,40,10,70,10,31,42])
        rdd = ranged_delta_decomposition(S,i=3,rv=(0.,100.),d=1)
        rdd1 = ranged_delta_decomposition(S,i=3,rv=(-20.,100.),d=-1)

        rdd_sol = [((0,3),-5),((3,4),-3),((4,21),-1),((21,30),1),((30,32),3),((32,60),5),((60,90),7)]
        rdd1_sol = [((0,-30),7)]

        assert rdd == rdd_sol
        assert rdd1 == rdd1_sol 

        S2 = np.array([13.5,14.9,39.1,10.15,69.7,8.1,33.4,44.4]) 
        rdd2 = ranged_delta_decomposition(S2,i=3,rv=(0.,100.),d=1)
        rdd3 = ranged_delta_decomposition(S2,i=3,rv=(-20.,100.),d=-1)

        rdd2_sol = [((0,3.35),-5),((3.35,4.75),-3),\
                ((4.75,23.25),-1),((23.25,28.95),1),\
                ((28.95,34.25),3),((34.25,59.55),5),\
                ((59.55,89.85),7)]
        assert rdd2 == rdd2_sol 

        rdd3_sol = [((0,-2.05),5),((-2.05,-30.15),7)] 
        assert rdd3 == rdd3_sol 

    def test__adjust_for_uwpd_change__case1(self): 
        S = np.array([13,14,40,10,70,10,31,42])
        S_ = adjust_for_uwpd_change(S,i=3,c=50,rv=(0.,100.),d_priority=1,recurse=True)
        S__ = adjust_for_uwpd_change(S,i=3,c=50,rv=(0.,100.),d_priority=-1,recurse=False)
        S1 = adjust_for_uwpd_change(S,i=3,c=100,rv=(0.,100.),d_priority=1,recurse=True)
        S2 = adjust_for_uwpd_change(S,i=3,c=-25,rv=(0.,100.),d_priority=1,recurse=True)

        pd0 = uwpd(S,accum_op=lambda x1,x2: x1 + x2)
        pd0_ = uwpd(S_,accum_op=lambda x1,x2: x1 + x2)
        assert pd0 + 50 == pd0_ 

        pd0__ = uwpd(S__,accum_op=lambda x1,x2: x1 + x2)
        assert pd0 + 50 == pd0__

        pd1 = uwpd(S1,accum_op=lambda x1,x2: x1 + x2)
        assert pd0 + 100 == pd1

        pd2 = uwpd(S2,accum_op=lambda x1,x2: x1 + x2)
        assert pd0 -25 == pd2
        
if __name__ == '__main__':
    unittest.main()