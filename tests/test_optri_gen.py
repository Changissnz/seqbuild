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

    def test__OpTriFlipDerivation__construct(self):
        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)
        otfd = OpTriFlipDerivation(ot,2,add,0)

        otfd.construct() 
        sol = np.array([[-21, -23, -24, -19, -17],\
            [  0, -29, -28, -22, -25],\
            [  0,   0, -19, -14, -23],\
            [  0,   0,   0, -12, -26],\
            [  0,   0,   0,   0, -10]], dtype=np.int32)
        assert np.all(otfd.m_ == sol) 

        otfd.reset_axis(1)
        otfd.construct() 
        sol = np.array([[ 10,   6,   8,   9,   4],\
            [  0,   8,  14,  13,   7],\
            [  0,   0,  11,   4,  -1],\
            [  0,   0,   0,  -5,  -3],\
            [  0,   0,   0,   0, -17]], dtype=np.int32)
        assert np.all(otfd.m_ == sol)


if __name__ == '__main__':
    unittest.main()