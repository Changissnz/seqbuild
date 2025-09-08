from seqgen.ssi_netop import * 
from morebs2.matrix_methods import is_proper_bounds_vector
import unittest

### lone file test 
"""
python -m tests.test_ssi_netop
"""
###
class SSINetOpMethods(unittest.TestCase):

    def test__SSINetOp__one_instance__case1(self):

        num_nodes = 20 
        prg = prg__LCG(56,93,176,901909)
        prg2 = prg__LCG(5,3,16,901)
        sno = SSINetOp.one_instance(num_nodes,prg,prg2)         
        assert True 

if __name__ == '__main__':
    unittest.main()