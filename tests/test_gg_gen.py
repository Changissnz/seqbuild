from seqgen.gg_gen import * 
from morebs2.numerical_generator import * 
import unittest


def GaugeGuidedGen__sample_XST(base_prg,target_measure): 
    base_output_span = [22,56] 
    density_measure = 45 
    super_range = [-100,2000]
    ngram_length = 5

    allow_subrange_drift = True   

    return GaugeGuidedGen(base_prg,base_output_span,\
        density_measure,super_range,ngram_length,target_measure,\
        allow_subrange_drift)



### lone file test 
"""
py -m tests.test_gg_gen
"""
###
class GaugeGuidedGenMethods(unittest.TestCase):

    def test__GaugeGuidedGen__next__case1(self):
        base_prg = prg__constant(3) 
        
        target_measure = "cov" 
        gg = GaugeGuidedGen__sample_XST(base_prg,target_measure)

        L = [next(gg) for _ in range(1000)] 
        l = len(set(L))
        assert l == 75, "got {}".format(l) 
        
        target_measure = "uwpd"
        gg = GaugeGuidedGen__sample_XST(base_prg,target_measure)

        L = [next(gg) for _ in range(1000)] 
        l = len(set(L))
        assert l == 50, "got {}".format(l) 

    def test__GaugeGuidedGen__next__case2(self): 
        #return 
        base_prg = prg__n_ary_alternator(-100,-50,-100)

        target_measure = "cov" 
        gg = GaugeGuidedGen__sample_XST(base_prg,target_measure)

        #L = [next(gg) for _ in range(1000)] 
        L = [] 
        for _ in range(1000): 
            x = next(gg) 
            print(x)
            L.append(x)
        l = len(set(L))
        assert l == 675, "got {}".format(l) 
        
        target_measure = "uwpd"
        gg = GaugeGuidedGen__sample_XST(base_prg,target_measure)

        L = [] 
        for _ in range(1000): 
            x = next(gg) 
            print(x)
            L.append(x)

        #L = [next(gg) for _ in range(1000)] 
        l = len(set(L))
        assert l == 255, "got {}".format(l) 


    def test__GaugeGuidedGen__next__case3(self): 
        #return 
        base_prg = prg__LCG(567,12,51,5001) 

        target_measure = "cov" 
        gg = GaugeGuidedGen__sample_XST(base_prg,target_measure)

        #L = [next(gg) for _ in range(1000)] 
        L = [] 
        for _ in range(1000): 
            x = next(gg) 
            print(x)
            L.append(x)
        l = len(set(L))
        assert l == 1000, "got {}".format(l) 
        
        target_measure = "uwpd"
        gg = GaugeGuidedGen__sample_XST(base_prg,target_measure)

        L = [] 
        for _ in range(1000): 
            x = next(gg) 
            print(x)
            L.append(x)

        l = len(set(L))
        assert l == 831, "got {}".format(l) 

    

if __name__ == '__main__':
    unittest.main()