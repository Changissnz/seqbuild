from desi.cy_metric import * 
from intigers.mod_prng import * 
from morebs2.numerical_generator import * 
import unittest

### lone file test 
"""
py -m tests.test_cy_metric 
"""
###
class CycleMetricMethods(unittest.TestCase):

    def test__CycleMetric__process_AND_summarize__case1(self):
 
        prg = prg__LCG(12,3,2,730) ## 54  

        num_iter = 150 
        pre_size = 3
        post_size = 3 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.11204, 20.13889])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.11204, 20.13889])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],5) == np.array([0.005, 0.5])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.11204, 20.13889])), "got {}".format(S['l-seq'])

        assert s0 == 1 
        assert s1 == 0.0 

        return 

    def test__CycleMetric__process_AND_summarize__case2(self):
 
        prg = prg__LCG(45.55,32.11,-345.776,3145.66) 
        prg = prg__single_to_int(prg) 

        num_iter = 150 
        pre_size = 3
        post_size = 3 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.009, 1.6])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.009, 1.6])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],5) == np.array([0.005, 0.5])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.009, 1.6])), "got {}".format(S['l-seq'])

        assert s0 == float('inf'), "got {}".format(s0)
        assert s1 == 0  
        return 

    def test__CycleMetric__process_AND_summarize__case3(self):
 
        prg = prg__iterable([4,5,6,7,8,9,10,12,13,14,15,16,17]) 

        num_iter = 150 
        pre_size = 3
        post_size = 3 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.08748,18.58974])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.08748, 18.58974])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],5) == np.array([0.005, 0.5])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.08748, 18.58974])), "got {}".format(S['l-seq'])

        assert s0 == 1., "got {}".format(s0)
        assert s1 == 0  
        return 

    def test__CycleMetric__process_AND_summarize__case4(self):
 
        prg = prg__iterable([4,5,6,7,8] + [10,4,6,5,12,15]) 

        num_iter = 150 
        pre_size = 1
        post_size = 1 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.01625, 4.875])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.01625, 4.875])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],5) == np.array([0.00799,1.22917])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.01625, 4.875])), "got {}".format(S['l-seq'])

        assert s0 == 1.375, "got {}".format(s0)
        assert s1 == 0.1875, "got {}".format(s1)
        return 

    def test__CycleMetric__process_AND_summarize__case5(self):
 
        prg = prg__iterable([3,13,76,54] + [54,76,13,13,3])
        
        num_iter = 150 
        pre_size = 1
        post_size = 1 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.11236, 18.83333])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.11236, 18.83333])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],5) == np.array([0.01481,2.72222])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.11236, 18.83333])), "got {}".format(S['l-seq'])

        assert s0 == 2.25, "got {}".format(s0)
        assert round(s1,5) == 0.36111, "got {}".format(s1)  
        
    def test__CycleMetric__process_AND_summarize__case6(self):
 
        prg = prg__LCG(525,562.4,-234,1003.5) 
        prg2 = prg__LCG(-5.25,156.24,-123.4,-703.5) 
        prg = merge_two_prgs(prg,prg2,add) 

        num_iter = 150 
        pre_size = 2
        post_size = 3 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.009, 1.6])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.009, 1.6])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],6) == np.array([0.005, 0.5])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.009, 1.6])), "got {}".format(S['l-seq'])

        assert s0 == float('inf'), "got {}".format(s0)
        assert s1 == 0  

        return 

    def test__CycleMetric__process_AND_summarize__case7(self):
 
        prg = prg__LCG(3,4,5,200) 
        prg2 = prg__LCG(7,5,-4,400) 
        prg = merge_two_prgs(prg,prg2,add) 

        num_iter = 150 
        pre_size = 2
        post_size = 3 

        M = CycleMetric(prg,num_iter,pre_size,post_size,cast_func=cr)
        M.process()
        S,s0,s1 = M.summarize() 

        assert np.all(np.round(S['seq'],5) == np.array([0.10576, 20.26667])), "got {}".format(S['seq'])
        assert np.all(np.round(S['set'],5) == np.array([0.10576, 20.26667])), "got {}".format(S['set'])
        assert np.all(np.round(S['valfreq'],6) == np.array([0.005, 0.5])), "got {}".format(S['valfreq'])
        assert np.all(np.round(S['l-seq'],5) == np.array([0.10576, 20.26667])), "got {}".format(S['l-seq'])

        assert s0 == float('inf'), "got {}".format(s0)
        assert s1 == 0  

        return 

if __name__ == '__main__':
    unittest.main()
