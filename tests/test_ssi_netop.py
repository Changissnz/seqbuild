from seqgen.ssi_netop import * 
from morebs2.matrix_methods import is_proper_bounds_vector
import unittest
import time 

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

    def test__SSINetOp__next__case1(self): 
        T0 = time.time() 

        num_nodes = 20 
        prg = prg__LCG(56,93,176,90199)
        prg2 = prg__LCG(5,3,16,901)

        # case 1 
        sno = SSINetOp.one_instance(num_nodes,deepcopy(prg),deepcopy(prg2),lcg_input_only=0,uniform_io_dist=1,shuffle_dist=0)
        sno.verbose = 0 

        sl = sno.struct_list

        Q0 = [] 
        for i in range(2000):
            q = next(sno)
            Q0.append(q)
        Q0 = np.array(Q0)

        # case 2 
        sno2 = SSINetOp.one_instance(num_nodes,deepcopy(prg),deepcopy(prg2),lcg_input_only=0,uniform_io_dist=0,shuffle_dist=1) 
        sno2.verbose = 0 

        Q1 = [] 
        for i in range(200): 
            q = next(sno2)
            Q1.append(q)
        Q1 = np.array(Q1)

        # case 3 
        sno3 = SSINetOp.one_instance(num_nodes,deepcopy(prg),deepcopy(prg2),lcg_input_only=1,uniform_io_dist=1,shuffle_dist=0) 
        sno3.verbose = 0 

        Q2 = [] 
        for i in range(2000): 
            q = next(sno3)
            Q2.append(q)
        Q2 = np.array(Q2)

        DX0 = Q0[:200] - Q1 
        i1 = np.where(DX0 != 0)[0]
        assert len(i1) == 196 

        DX1 = Q2 - Q0 
        i2 = np.where(DX1 != 0)[0]
        assert len(i2) == 50,"got {}".format(len(i2))

        T1 = time.time() - T0 
        print("elapsed time: {}".format(T1))

if __name__ == '__main__':
    unittest.main()