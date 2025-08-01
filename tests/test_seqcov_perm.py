from desi.seqcov_perm import *
from morebs2.numerical_generator import prg__LCG
import unittest

### lone file test 
"""
python -m tests.test_seqcov_perm
"""
###
class SeqCoveragePermuterMethods(unittest.TestCase):

    '''
    test demonstrates sensitivity to changes in input, 
    via comparison of the two "very" different output vectors, 
    despite the only difference in index-to-index comparison of 
    input vectors being index 4. These changes could not have been 
    possible under linear vector operations (normal,index-to-index), 
    such as multiplication. 
    '''
    def test__SeqCoveragePermuter__apply_posANDneg_delta__case1(self):

        sequence = np.array([-3,-1,0.5,0.1,7,-10,10,12,13,17])
        coverage_delta = 0.2 
        max_radius = 0.5 
        super_range = (-10,20.)
        prg = prg__LCG(78,112,-13,1112)
        scp = SeqCoveragePermuter(sequence,coverage_delta,max_radius,super_range,prg)
        scp.set_partition()

        cov00 = scp.cov_typeabs

        scp.apply_pos_delta()
        cov01 = scp.cov_typeabs

        assert cov00 == 8.9
        assert cov01 == 14.9

        rs = np.array([\
        [-10.,-9.5],\
        [-3.5,-2.5],\
        [-1.5,-0.5],\
        [-0.4,1.3],\
        [6.5,7.5],\
        [9.5,10.5],\
        [11.5,13.5],\
        [ 16.5,17.5]]) 

        crs = np.array([\
            [-9.5, -3.5],\
            [-2.5, -1.5],\
            [ 1.3, -0.4],\
            [ 1.,   6.5],\
            [ 7.5,  9.5],\
            [10.5, 11.5],\
            [13.5, 16.5],\
            [17.5, 20. ]]) 

        changean = 1.6 
        coverage_delta2 = coverage_delta * -1 
        scp2 = SeqCoveragePermuter(sequence,coverage_delta2,max_radius,super_range,prg)
        scp2.set_partition()
        cov10 = scp2.cov_typeabs

        qx = deepcopy(scp2.rs) 
        scp2.apply_neg_delta()
        cov11 = scp2.cov_typeabs

        assert cov11 + 4.0 == cov10 

if __name__ == '__main__':
    unittest.main()



