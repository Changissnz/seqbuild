from intigers.process_seq import * 
import unittest

### lone file test 
"""
python3 -m tests.test_seq2graph
"""
###
class ProcessSeqMethods(unittest.TestCase):

    def test__gleqvec(self):
        vec1 = np.array([2,4,8,16,3,6,12,3,9,27,81],dtype=np.int32)
        vec2 = np.array([2,0,4,8,16,0,3,6,12,0,3,9,27,81],dtype=np.int32)
        vec3 = np.array([0,0,2,2,2.5514132,2.55142],dtype=np.float32)

        gvec = gleqvec(vec1,rounding_depth=5)
        sol_gvec = np.array([ 1,  1,  1, -1,  1,  1, -1,  1,  1,  1], dtype=np.int32)
        assert (gvec == sol_gvec).all() 

        gvec2 = gleqvec(vec2,rounding_depth=5)
        sol_gvec2 = np.array([-1,  1,  1,  1, -1,  1,  1,  1, -1,  1,  1,  1,  1], dtype=np.int32)
        assert (gvec2 == sol_gvec2).all()

        gvec3 = gleqvec(vec3,rounding_depth=5)
        sol_gvec3 = np.array([0, 1, 0, 1, 1], dtype=np.int32)
        assert (gvec3 == sol_gvec3).all() 


if __name__ == '__main__':
    unittest.main()