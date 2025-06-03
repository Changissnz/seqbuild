from seqgen.rch_gen import * 
from morebs2.numerical_generator import prg__LCG 

import unittest

### lone file test 
"""
python -m tests.test_rch_gen
"""
###
class RCHGenMethods(unittest.TestCase):

    def test__RCHAccuGen__one_new_RCHAccuGen__v1(self):
        num_nodes = 2
        dim_range = [3,6]
        prg = prg__LCG(3,4,5,33)
        ufreq_range = [2,3]
        mutrate = 1.0 
        queue_capacity = 1000 

        rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
                ufreq_range,mutrate,queue_capacity)

        # size check with |<RChainHead>|
        rg.apply(16) 
        assert len(rg.acc_queue) == 3 

        rg.apply(30)
        rg.apply(233) 
        rg.apply(52) 

        # counter check for update 
        D = {0: {'cf': 2, 'rf': 2}, 1: {'rf': 2, 'cf': 2}}
        assert rg.update_log == D 

        # size check for queue
        for _ in range(1200): 
            rg.apply(prg())
        assert len(rg.acc_queue) == 1000 

if __name__ == '__main__':
    unittest.main()