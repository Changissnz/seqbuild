from seqgen.udls_gen import * 
from morebs2.numerical_generator import * 
import unittest

### lone file test 
"""
py -m tests.test_udls_gen
"""
###
class UDLSGenMethods(unittest.TestCase):

    def test__UDLSGen__next__case1(self): 

        prg = prg__LCG(56.66,-17.7,321.56,10103.55) 
        varsize_range = [4,8]
        allow_rate1_change = True 
        allow_rate2_change = False 

        ug = UDLSGen(prg,varsize_range,allow_rate1_change,allow_rate2_change,\
            coeff_absmax=DEFAULT_UDLSGEN_COEFF_ABSMAX,\
            draw_values_from_cache=True)

        L = [next(ug) for _ in range(5000)]
        S = set(np.round(L,5)) 
        assert len(S) == 5000 

    def test__UDLSGen__next__case2(self):

        prg = prg__LCG(56.66,-17.7,321.56,10103.55) 
        varsize_range = [4,8]
        allow_rate1_change = True 
        allow_rate2_change = False 

        ug = UDLSGen(deepcopy(prg),varsize_range,allow_rate1_change,allow_rate2_change,\
            coeff_absmax=DEFAULT_UDLSGEN_COEFF_ABSMAX,\
            draw_values_from_cache=True)
        ug2 = UDLSGen(deepcopy(prg),varsize_range,allow_rate1_change,allow_rate2_change,\
            coeff_absmax=DEFAULT_UDLSGEN_COEFF_ABSMAX,\
            draw_values_from_cache=False)
        ug3 = UDLSGen(deepcopy(prg),varsize_range,allow_rate1_change,allow_rate2_change=True,\
            coeff_absmax=DEFAULT_UDLSGEN_COEFF_ABSMAX,\
            draw_values_from_cache=True)

        L = [next(ug) for _ in range(5000)]
        L2 = [next(ug2) for _ in range(5000)]
        L3 = [next(ug3) for _ in range(5000)]        

        S = set(np.round(L,5)) 
        S2 = set(np.round(L2,5)) 
        S3 = set(np.round(L3,5)) 

        assert len(S2.intersection(S)) == 182
        assert len(S3.intersection(S)) == 0, "got {}".format(len(S3.intersection(S)))
        assert len(S3.intersection(S2)) == 0, "got {}".format(len(S3.intersection(S)))


if __name__ == '__main__':
    unittest.main()
