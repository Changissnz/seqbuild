from mini_dm.sp_op import * 
from morebs2.numerical_generator import prg__LCG
from mini_dm.minmax_freq import * 
import unittest

### lone file test 
"""
python -m tests.test_sp_op
"""
###
class SuperPartitionOpMethods(unittest.TestCase):

    def test__SuperPartitionOp__one_partition__case1(self):

        prg = prg__LCG(896,-452,60041,100523)
        sp2 = [set([0,1,2,3,4,5,6,7,8]),\
            set([1,9,10]),\
                set([2,9,11]),\
                set([10,11]),\
                set([0,1,2,3]),\
                set([4,5,6,7,8]),\
                set([12,13,14]),\
                set([12]),\
                set([13,14,15,16,17]),\
                set([16,17,18,19,20]),\
                set([18,19,20])] 
        spo2 = SuperPartitionOp(sp2,prg)  

        qx = []
        for _ in range(5): 
            p = spo2.one_partition() 
            if p in qx: 
                continue 
            qx.append(p) 
        assert len(qx) == 5 

        rx = flatten_setseq(sp2) 
        for qx_ in qx: 
            assert flatten_setseq(qx_) == rx 
            fm = setseq_to_frequency_map(qx_)
            vfm = set(fm.values())
            assert vfm == {1}

        return 
    

if __name__ == '__main__':
    unittest.main()