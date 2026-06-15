from seqgen.ssi_netop import * 
from morebs2.matrix_methods import is_proper_bounds_vector
import unittest
import time 

"""
An <SSINetOp> with 0 LCG nodes. 
"""
def SSINetOp__sample_0LCG(rapid_update:bool): 

    prg = prg__LCG(56,93,176,9019)
    batch_size = 15 

    sl0 = SSINetOp.one_instance__slist("mdr",batch_size,prg)
    sl1 = SSINetOp.one_instance__slist("optri",batch_size,prg)

    sl0.extend(sl1) 

    num_trees = 3 
    ssb2n = SSIBatch2Net(sl0,num_trees,prg) 
    ssb2n.make_net()
    h2t = ssb2n.h2tree_map

    return SSINetOp(sl0,h2t,prg,lcg_input_only=0,uniform_io_dist=1,\
        shuffle_dist=0,io_prg=None,rapid_update=rapid_update,verbose=False) 


### lone file test 
"""
py -m tests.test_ssi_netop
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
        assert len(i1) == 194, "got {}".format(len(i1)) 

        DX1 = Q2 - Q0 
        i2 = np.where(DX1 != 0)[0]
        assert len(i2) == 1968,"got {}".format(len(i2))

        T1 = time.time() - T0 
        print("elapsed time: {}".format(T1))

    """
    compares difference b/t rapid and slow updating of SSINet node 
    parameters
    """ 
    def test__SSINetOp__next__case2(self): 

        rapid_ssn = SSINetOp__sample_0LCG(rapid_update= True) 
        slow_ssn = SSINetOp__sample_0LCG(rapid_update= False) 

        # run the rapid-updating 
        t = time.time() 
        for _ in range(25000): 
            next(rapid_ssn) 
            if _ % 250 == 0: 
                print("iter: ",_)

        t0 = time.time() - t 

        M0 = rapid_ssn.node_to_activation_count_map() 
        assert M0 == {0: 1230, 1: 1192, 2: 548, 3: 597, 4: 1, \
            5: 1531, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, \
            11: 1, 12: 1, 13: 1, 14: 1, 15: 1107, \
            16: 507, 17: 1, 18: 1, 19: 1, 20: 1, \
            21: 1, 22: 1, 23: 1, 24: 1, 25: 1, \
            26: 1, 27: 1, 28: 1, 29: 1}

        # run the slow-updating 
        t = time.time() 
        for _ in range(25000): 
            next(slow_ssn) 
            if _ % 250 == 0: 
                print("iter: ",_)

        t1 = time.time() - t 

        M1 = slow_ssn.node_to_activation_count_map() 
        assert M1 == {0: 38, 1: 26, 2: 15, 3: 28, 4: 1, \
            5: 43, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, \
            11: 1, 12: 1, 13: 1, 14: 1, 15: 3, \
            16: 5, 17: 1, 18: 1, 19: 1, 20: 1, \
            21: 1, 22: 1, 23: 1, 24: 1, 25: 1, \
            26: 1, 27: 1, 28: 1, 29: 1}

        # check that slow-updating is faster 
        assert t1 < t0 
        return

if __name__ == '__main__':
    unittest.main()