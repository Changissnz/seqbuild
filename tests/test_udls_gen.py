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


if __name__ == '__main__':
    unittest.main()
