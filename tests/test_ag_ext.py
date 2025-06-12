from intigers.ag_ext import * 
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
        
if __name__ == '__main__':
    unittest.main()