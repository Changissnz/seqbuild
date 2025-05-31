from intigers.idt_proc import * 
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
        xr = is2t.split__depthreq(tn,"factor")

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

    def test__IntSeq2Tree__poly_split__depthreq__case1(self): 
            
        prng = prg__constant(0)
        l = None 
        d = 4 
        L = [3,4,10,20,32,12,17,15,16,2,1,5,13]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng)
        is2t.init_root() 
        tn = is2t.node_cache.pop(0) 
        q = is2t.split__depthreq(tn,"poly")
        tn = is2t.node_cache[0] 

        D = defaultdict(int)
        for l in L: 
            q = tn.travf.apply(l)
            D[q] += 1 
        assert D == {None:9,1:4}

        c = 0 
        while True: 
            if len(tn.children) == 0: 
                break     
            tn = tn.children[0]
            c += 1 
        assert c == 4

    def test__IntSeq2Tree__poly_subset_bclassifier(self): 
        prng = prg__constant(0)
        l = None 
        d = 4 
        L = [480,320,6400,1280,804,7,19,11,50000,735,230,170,556]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng)

        q,_ = is2t.poly_subset_bclassifier(deepcopy(L),4)

        r = 0 
        for l in L: 
            if q.bclassify(l): r += 1 
        assert r == 4 

    """
    simple case of integer sequence, relatively small numbers, 
    no duplicates, no zeros. 
    """
    def test__IntSeq2Tree__convert__case1(self): 
        prng = prg__constant(0)#3)
        l = None 
        d = 4 
        L = [123, 321, 43, 15, 17, 18, 19]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng,verbose=False)
        is2t.convert()
        q = is2t.root
        dr,depth = TNode.dfs(q,display=False,collect=True,reset_index=True)

        D = {}
        for x in L: 
            y = q.travf.apply(x)
            D[x] = y 

        assert depth == 4
        assert D == {123: 1, 321: 1, \
            43: 1, 15: 1, 17: 5, 18: 5, 19: 6}

        itp = IDTProc(q)
        D2 = dict() 
        D2[123] = [0, 1]
        D2[321] = [0, 1, 2]
        D2[43] = [0, 1, 2, 3]
        D2[15] = [0, 1, 2, 3, 4]
        D2[17] = [0, 5, 8]
        D2[18] = [0, 5, 7]
        D2[19] = [0, 6]

        for l in L[:7]: 
            p = itp.process_value(l)
            assert D2[l] == p 

if __name__ == '__main__':
    unittest.main()