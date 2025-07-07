from mini_dm.sp_op import * 
from morebs2.numerical_generator import prg__LCG
from mini_dm.minmax_freq import * 
from mini_dm.vs_fit_op import * 
from morebs2.search_space_iterator import * 
import numpy as np 
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
            p,_ = spo2.one_partition() 
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

    """
    same example as that of method<test__IOFit__superpart_by_AffineDelta__case1>. 
    """
    def test__SuperPartitionOp__smallestANDlargest_partition(self):
        m = np.array([4,5,16,27,54,65]) 
        a = 43.0 
        ad = AffineDelta(m,a,0) 

        bounds = np.array([[-1.0,1.0],\
                        [2.0,7.0],\
                            [14.0,22.0],\
                            [-20,-2],\
                            [55.,105.0],\
                            [200.,490.]])

        start_point = deepcopy(bounds[:,0])
        column_order = np.arange(6)
        ssi_hop = np.array([2,5,8,18,50,290]) 
        cycle_on = True 
        cycle_is = 0 

        ssi = SearchSpaceIterator(bounds,start_point,\
            column_order,ssi_hop,cycle_on,cycle_is)
        xs,ys = [],[]

        for i in range(100):
            x = next(ssi) 
            xs.append(x) 
            ys.append(ad.fit(x))

        iof = IOFit(xs,ys,None,None,None)
        sprt = iof.superpart_by_AffineDelta(0)

        prg = prg__LCG(896,-452,60041,100523)
        sp = deepcopy(sprt.partition)
        spo = SuperPartitionOp(sp,prg)

        qx = spo.smallest_partition() 
        qx2 = spo.largest_partition() 
        assert len(qx[0]) == 1 
        assert len(qx2[0]) == 57 


if __name__ == '__main__':
    unittest.main()