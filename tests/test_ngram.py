from mini_dm.ngram import * 
from morebs2.numerical_generator import prg__LCG
import unittest

### lone file test 
"""
python3 -m tests.test_ngram 
"""
###
class NGrammer2DMethods(unittest.TestCase):

    """
    case demonstrates improvement of replication 
    from <ModuloDecomp> to <ModuloDecompV2>. 
    """
    def test__NGrammer2D__one_cycle__case1(self):
        L = [10,20,30,40,50,60,70,8,9,10,11] 
        n = 4
        dim2 = (2,2) 
        ng = NGrammer2D(L,n,dim2)
        c = list(ng.one_cycle()) 
        assert len(c) == 11

        sol0 = np.array([\
            [10, 20],\
            [30, 40]]),

        sol1 = np.array([\
            [11, 10],\
            [20, 30]])

        assert np.all(c[0] == sol0) 
        assert np.all(c[-1] == sol1) 
        return


if __name__ == '__main__':
    unittest.main()


