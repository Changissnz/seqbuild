from intigers.idt_proc import * 
from morebs2.numerical_generator import prg__constant,\
    prg__n_ary_alternator,LCG,prg__LCG
from morebs2.globalls import std_invert_map 

import unittest

### lone file test 
"""
python -m tests.test_idectree 
"""
###
class IDecTreeMethods(unittest.TestCase):

    def test__partition_fix__subset_is_minsize_2__case1(self):
        g = [1,4,1] 
        q = prg__n_ary_alternator(-10,15,5)
        g2 = partition_fix__subset_is_minsize_2(g,q)
        assert np.all(g2[0] == [2,2,2])
        assert g2[1] == True

        g = [1,3,1] 
        g2 = partition_fix__subset_is_minsize_2(g,q)
        assert g2[1] == False

        g = [30,1,1,1,1,1,1,27,1,1,1,1,1] 
        g2 = partition_fix__subset_is_minsize_2(g,q)
        assert np.all(g2[0] == [24,  2,  2,  2,  2,  \
            2,  2, 22,  2,  2,  2,  2,  2])
        assert g2[1] == True

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
        assert dx == {None:1,2:3}, "got {}".format(dx) 

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

    Tests for depth requirement, correct labelling by `travf`. 
    """
    def test__IntSeq2Tree__convert__case1(self): 
        prng = prg__constant(0)
        l = None 
        d = 4 
        L = [123, 321, 43, 15, 17, 18, 19]

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng,verbose=False)
        is2t.convert()
        q = is2t.root
        dr,depth,_ = TNode.dfs(q,display=False,collect=True,reset_index=True)

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

        for l in L: 
            p,_ = itp.process_value(l)
            assert D2[l] == p 

    """
    Similar test case to case 1. Uses larger values than case 1 and 
    demonstrates with two constant generators.
    """
    def test__IntSeq2Tree__convert__case2(self): 
        l = None 
        d = 4 

        # prng #1 
        prng = prg__constant(0)#3)

        L = [480,320,6400,1280,804,7,19,482,330,350,6000,1220,1032,450] 
        is2t = IntSeq2Tree(IntSeq(L),l,d,prng,verbose=False)
        is2t.convert()

        q = is2t.root
        dr,depth,_ = TNode.dfs(q,display=False,collect=True,reset_index=True)
        assert depth == 4 

        D = {480: 1, 320: 1, 6400: 1, 1280: 1, 804: 5, \
            7: 5, 19: 6, 482: 5, 330: 6, 350: 6, 6000: 6, \
            1220: 5, 1032: 6, 450: 5}

        for x in L:
            y = q.travf.apply(x)
            assert D[x] == y 

        D2 = std_invert_map(D) 
        assert len(D2[1]) == 4   

        # prng #2 
        prng = prg__constant(5)#3)
        is2t = IntSeq2Tree(IntSeq(L),l,d,prng,verbose=False)
        is2t.convert()

        q = is2t.root
        dr,depth,_ = TNode.dfs(q,display=False,collect=True,reset_index=True)
        assert depth == 4 

    """
    test for correct leaf size
    """
    def test__IntSeq2Tree__convert__case3(self): 
        prng = prg__n_ary_alternator(-8,8,-7)

        L = [48,32,640,128,804,7,19,482,330,350,600,1220,1032,450] 
        l = 5
        d = None 

        is2t = IntSeq2Tree(IntSeq(L),l,d,prng,verbose=False)
        is2t.convert()
        q = is2t.root

        assert len(q.children) == 5 

        prng = prg__LCG(-4,5,2,29)
        is2t = IntSeq2Tree(IntSeq(L),l,d,prng,verbose=False)
        is2t.convert()
        q = is2t.root
        assert len(q.children) == 5 

    """
    capacity test for input magnitude 
    """
    def test__IntSeq2Tree__convert__case4(self):     
        prg = prg__LCG(3,7,9,33)
        prg2 = prg__n_ary_alternator(-20,30,-19) 

        q = [] 
        for i in range(50): 
            q_ = prg()
            q2_ = prg2() 
            q.append((q_ * q2_) % 422)
        q = list(set(q) - {0})
        
        ns = IntSeq(q) 
        l = None#4 
        d = 6 
        is2t = IntSeq2Tree(ns,l,d,prg,verbose=True) 
        is2t.convert()

        ns = IntSeq(q) 
        l = 4 
        d = None  
        is2t = IntSeq2Tree(ns,l,d,prg,verbose=True) 
        is2t.convert()

    """
    capacity test #2 for input magnitude 
    """
    def test__IntSeq2Tree__convert__case5(self):
        from random import randrange 

        qx = np.unique([randrange(2,20000) for _ in range(100)])
        iq = IntSeq(qx)
        prg1 = prg__LCG(13,87,2011,500)
        is2t = IntSeq2Tree(iq,2,None,prg1,verbose=True) 
        is2t.convert()


if __name__ == '__main__':
    unittest.main()