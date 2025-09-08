from seqgen.ssi_load import * 
import unittest

### lone file test 
"""
python -m tests.test_ssi_load
"""
###

class OpTriGenLiteMethods(unittest.TestCase):

    def test__OpTriGenLite__init__case1(self):
        L = [i for i in range(4,12)] 
        L = IntSeq(L) 

        otl = OpTriGenLite(L,add,sub)
        otl_sol = np.array([\
            [ 1,  1,  1,  1,  1,  1,  1],
            [64,  0,  0,  0,  0,  0,  0],
            [32, 32,  0,  0,  0,  0,  0],
            [16, 16, 16,  0,  0,  0,  0],
            [ 8,  8,  8,  8,  0,  0,  0],
            [ 4,  4,  4,  4,  4,  0,  0],
            [ 2,  2,  2,  2,  2,  2,  0]]) 
        assert np.all(otl.m == otl_sol) 

class SSIBatchLoader__TypeLCGNetMethods(unittest.TestCase):

    def test__SSIBatchLoader__TypeLCGNet__instantiate_slist__case1(self):
        sidn = "lcg"
        param_bounds = np.array([\
            [-4,5],\
            [-4,5],\
            [-4,5],\
            [-4,-3.995]])

        max_batch_size = float('inf')
        aux_prg = prg__LCG(63,3,199,4000.0)

        ssibl = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,max_batch_size,aux_prg)

        ssibl.instantiate_slist()
        assert len(ssibl.slist) == 648

    def test__SSIBatchLoader__TypeLCGNet__instantiate_slist__case2(self):
        sidn = "mdr"
        param_bounds = (0,0)
        max_batch_size = 10

        aux_prg = prg__LCG(13,43,651,4120.0)
        ssibl = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,max_batch_size,aux_prg)

        ssibl.instantiate_slist()
        assert len(ssibl.slist) == 10 

    def test__SSIBatchLoader__TypeLCGNet__instantiate_slist__case3(self):
        sidn = "optri"
        param_bounds = None 
        max_batch_size = 7 

        aux_prg = prg__LCG(13,43,651,4120.0)
        ssibl = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,max_batch_size,aux_prg)

        ssibl.instantiate_slist()
        assert len(ssibl.slist) == 7 

class SSIBatch2NetMethods(unittest.TestCase): 

    def test__SSIBatch2Net__make_net__case1(self): 
        ## 
        sidn = "optri"
        param_bounds = None 
        max_batch_size = 4 

        aux_prg = prg__LCG(13,43,651,4120.0)
        ssibl = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,max_batch_size,aux_prg)

        ssibl.instantiate_slist()

        ## 
        sidn = "lcg"
        param_bounds = np.array([\
            [-4,5],\
            [-4,50],\
            [-4,5],\
            [-4,6]])

        max_batch_size = 4 
        aux_prg = prg__LCG(63,3,199,4000.0)

        ssibl2 = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,max_batch_size,aux_prg)
        ssibl2.instantiate_slist()

        ## 
        sidn = "mdr"
        param_bounds = (0,0)
        max_batch_size = 4

        aux_prg = prg__LCG(13,43,651,4120.0)
        ssibl3 = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,max_batch_size,aux_prg)

        ssibl3.instantiate_slist()
        slist = ssibl.slist + ssibl2.slist + ssibl3.slist 

        ssb2n = SSIBatch2Net(slist,5,aux_prg) 
        ssb2n.make_net()

        all_nodes = set([i for i in range(12)])
        for k in ssb2n.h2tree_map.keys():
            T = ssb2n.h2tree_map[k] 

            kx = set(T.keys()) 
            for v in T.values(): 
                kx = kx | v 
            assert kx == all_nodes

if __name__ == '__main__':
    unittest.main()