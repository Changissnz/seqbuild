from seqgen.optri_gen import * 
import unittest


### lone file test 
"""
python3 -m tests.test_optri_gen
"""
###

class OpTriGenMethods(unittest.TestCase):

    def test__OpTri45N90Split__j_derivative_seq(self):
        split = ([12,13,14],[15,16,17]) 
        ots = OpTri45N90Split(split)
        s = ots.j_derivative_seq(1)

        assert s == [12, 25, 52, 94, 152, 227]
        assert ots.derivative_seqs == \
            {3: [14, 15, 16, 17],\
            2: [13, 27, 42, 58, 75],\
            1: [12, 25, 52, 94, 152, 227]}


if __name__ == '__main__':
    unittest.main()