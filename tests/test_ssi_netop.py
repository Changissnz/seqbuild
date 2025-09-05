from seqgen.ssi_load import * 
from morebs2.matrix_methods import is_proper_bounds_vector
import unittest

### lone file test 
"""
python -m tests.test_ssi_netop
"""
###
class SSINetOpMethods(unittest.TestCase):

    def test__SSINetOp__TypeLCGNet__instantiate_slist__case1(self):
        aux_prg = prg__LCG(13,43,651,4120.0)
        aux_prg = wrap_ranged_modulo_over_generator(aux_prg,[-2050.0,2050.0])
        b_out = prg__single_to_bounds_outputter(aux_prg,6) 

        for i in range(5): 
            bd = b_out()
            assert is_proper_bounds_vector(bd) 

if __name__ == '__main__':
    unittest.main()