from intigers.idectree import * 
from morebs2.numerical_generator import prg__constant,prg__n_ary_alternator,LCG 

import unittest

### lone file test 
"""
python -m tests.test_idectree 
"""
###
class IDecTreeMethods(unittest.TestCase):

    def test__IntSeq2Tree__factor_split__partitioned__case1(self):
        prng = prg__constant(0)

        l = None 
        d = 4 
        L = [480,320,6400,1280,804,7,19,11]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng)
        SX = deepcopy(L) 
        partition = [4,2,2] 

        is2t2 = deepcopy(is2t) 
        q = is2t2.factor_split__partitioned(L,partition,last_subset_isneg=True)
        q2 = is2t.factor_split__partitioned(L,partition,last_subset_isneg=False) 

        assert q == {0: [(5, {480,6400,320,1280})], \
            1: [(134, {804}), (19, {19})]}

        assert q2 == {0: [(5, {480,6400,320,1280})], \
            1: [(134, {804}), (19, {19})], \
            2: [(1, {480, 6400, 320, 1280, 804, 7, 11, 19})]},"got {}".format(q2)
        
        #{0: [(5, {480,6400,320,1280})], \
        #    1: [(134, {804}), (19, {19})], 2: [(1, {11,7})]}, "got {}".format(q2)

    def test__IntSeq2Tree__factor_split__partitioned__case2(self):

        prng = prg__constant(0)

        l = None 
        d = 4 
        L = [480,320,6400,1280,804,7,19,11,50000,735,230,170,556]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng)
        SX = deepcopy(L) 
        partition = [4,2,2,len(L) - 8] 

        q2 = is2t.factor_split__partitioned(L,partition,last_subset_isneg=False) 

        assert q2 ==  {0: [(32, {480,6400,320,1280})], \
            1: [(7, {7,735})], 2: [(11, {11}), (268, {804})],\
             3: [(1, {480, 6400, 320, 1280, 804, 230, \
                7, 170, 11, 556, 50000, 19, 735})]}


if __name__ == '__main__':
    unittest.main()