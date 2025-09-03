from seqgen.ssi_load import * 
import unittest

### lone file test 
"""
python -m tests.test_ssi_load
"""
###
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

if __name__ == '__main__':
    unittest.main()