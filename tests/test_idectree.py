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
            2: [(1, {11, 7})]},"got {}".format(q2)

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
            1: [(7, {7,735})], 2: [(11, {11}), (268, {804})], \
            3: [(1, {230, 170, 556, 50000, 19})]}, "got {}".format(q2) 

    def test__IntSeq2Tree__factor_split__depthreq__case1(self): 

        prng = prg__constant(0)
        l = None 
        d = 4 
        L = [480,320,6400,1280,804,7,19,11,50000,735,230,170,556]
        is2t = IntSeq2Tree(IntSeq(L),l,d,prng)

        is2t.init_root() 
        tn = is2t.node_cache.pop(0) 
        xr = is2t.factor_split__depthreq(tn)

        assert len(tn.acc_queue) == 13 - 4
        for x in tn.acc_queue: 
            assert type(tn.travf.apply(x)) == type(None)


        tn2 = tn.children[0]
        elem = tn2.acc_queue
        dx = defaultdict(int)
        for e in elem:
            dx[tn2.travf.apply(e)] += 1 
        assert dx == {None:1,2:3} 

        tn3 = tn2.children[0]
        dx = defaultdict(int)
        for e in tn3.acc_queue:
            dx[tn3.travf.apply(e)] += 1 
        assert dx == {None:1,3:2} 

        tn4 = tn3.children[0]
        dx = defaultdict(int)
        for e in tn4.acc_queue:
            dx[tn4.travf.apply(e)] += 1 
        assert dx == {None:1,4:1}  

        tn5 = tn4.children[0]
        dx = defaultdict(int)
        assert type(tn5.travf) == type(None)
        return 

    def test__IntSeq2Tree__poly_subset_classifier(self): 
        prng = prg__constant(0)
        l = None 
        d = 4 
        L = [480,320,6400,1280,804,7,19,11,50000,735,230,170,556]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng)

        q = is2t.poly_subset_classifier(deepcopy(L),4)

        r = 0 
        for l in L: 
            if q.bclassify(l): r += 1 
        assert r == 4 

if __name__ == '__main__':
    unittest.main()