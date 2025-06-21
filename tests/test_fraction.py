from desi.fraction import * 
from morebs2.aprng_gauge import BatchIncrStruct
from morebs2.numerical_generator import prg__LCG,prg__constant
from intigers.extraneous import prg__LCG_to_range_outputter,prg__LCG_to_ndim_index_outputter
import unittest

def fraction__sample_LCGs(l):

    lx = prg__LCG(87,14,53,799) 
    rx = prg__LCG(166.341,-13.4567,5321.756267,5500.1543)
    ix = prg__LCG(25,3,87,2024)

    rxi = prg__LCG_to_range_outputter(rx) 
    outi = prg__LCG_to_ndim_index_outputter(ix,(l,l))
    return outi,lx,rxi 


def fraction__sample_LCGs_2(l): 
    lx = prg__constant(10)

    bis = BatchIncrStruct(l,is_perm=True,\
        is_reflective=True,subset_size=2)
    ix = bis.__next__ 

    rx = prg__constant((0.,5000.0))
    return ix,lx,rx 

### lone file test 
"""
python -m tests.test_fraction 
"""
###
class FractionMethods(unittest.TestCase):

    def test__QValueOutputter__next__case1(self):
        intseq0 = IntSeq([1,0,1,0,10,0,1,0,10,1,0])
        intseq = IntSeq([4567,234,15564,1134,51,78,95547,46892,1416,747345223])

        outi,lx,rxi = fraction__sample_LCGs(len(intseq)) 
        qvo = QValueOutputter(intseq,outi,lx,rxi,1) 
        qvo2 = QValueOutputter(intseq0,outi,lx,rxi,1) 

        l0,l2 = [],[]
        for i in range(5): 
            l0.append(next(qvo))
            l2.append(next(qvo2)) 

        q1 = [171.341, 2367.146623956911, 5632.743898264879, 4608.298404089755, 4005.4265540921688]
        q2 = [2331.0486668385947, 4450.995057997054, 1752.823161311726, 338.44678641201244, 4425.551978933781]

        # value test
        assert np.all(q1 == l0)
        assert np.all(q2 == l2)

        outi,lx,rxi = fraction__sample_LCGs(len(intseq)) 
        qvo_ = QValueOutputter(intseq,outi,lx,rxi,1) 
        qvo2_ = QValueOutputter(intseq0,outi,lx,rxi,1) 

        l1 = []
        for i in range(5):
            l1.append(next(qvo_))
            next(qvo2_)

        # determinism test 
        assert np.all(l1 == q1)

        qvo_.change_adj_type()
        l3,l4 = [],[] 
        for i in range(5): 
            l3.append(next(qvo_)) 
            l4.append(next(qvo)) 

        # different adjustment type, different output test 
        assert not np.any(l3 == l4) 

    def test__QValueOutputter__next__case2(self):

        intseq = IntSeq([1,2,4,2,1,4,1,2,2,4])
        outi,lx,rxi = fraction__sample_LCGs_2(len(intseq))

        qvo = QValueOutputter(intseq,outi,lx,rxi,1) 
        l = [] 
        for i in range(10): 
            l.append(next(qvo))
        assert np.all(l == [0.0, 5.0, 25.0, 5.0, 0.0, 25.0, 0.0, 5.0, 5.0, 25.0]) 


if __name__ == '__main__':
    unittest.main()