from intigers.seq_struct import * 
import unittest

### lone file test 
"""
python3 -m tests.test_seq_struct
"""
###
class SeqStructMethods(unittest.TestCase):

    def test__IntSeq__difftri(self):

        its = IntSeq([5,4,3])
        prod = its.difftri() 
        assert (prod == np.array([[-1, -1],[ 0,  0]], dtype='int32')).all()

        its2 = IntSeq([15,40,312])
        prod2 = its2.difftri() 
        assert (prod2 == np.array([[ 25, 272],[  0, 247]], dtype='int32')).all()

        its3 = IntSeq([2,4,8,16,3,6,12,3,9,27,81])
        prod3 = its3.difftri() 
        sol3 = np.array([[   2,    4,    8,  -13,    3,    6,   -9,    6,   18,   54],\
            [   0,    2,    4,  -21,   16,    3,  -15,   15,   12,   36],\
            [   0,    0,    2,  -25,   37,  -13,  -18,   30,   -3,   24],\
            [   0,    0,    0,  -27,   62,  -50,   -5,   48,  -33,   27],\
            [   0,    0,    0,    0,   89, -112,   45,   53,  -81,   60],\
            [   0,    0,    0,    0,    0, -201,  157,    8, -134,  141],\
            [   0,    0,    0,    0,    0,    0,  358, -149, -142,  275],\
            [   0,    0,    0,    0,    0,    0,    0, -507,    7,  417],\
            [   0,    0,    0,    0,    0,    0,    0,    0,  514,  410],\
            [   0,    0,    0,    0,    0,    0,    0,    0,    0, -104]],\
            dtype=np.int32)
        assert (prod3 == sol3).all()

    def test__AffineFitSearch__count(self): 
        l = [2,4,8,16,3,6,12,3,9,27,81]
        afs = AffineFitSearch(l,exclude_neg=True,log_revd=True)

        afs.load_all_candidates()
        afs.count()
        assert afs.mmf.sorted_counts[-3:] == [((1, 6), 2), ((3, 0), 3), ((2, 0), 5)]

    def test__AffineFitSearch__default_affine_decomp(self):
        l = [2,4,8,16,3,6,12,3,9,27,81]
        afs = AffineFitSearch(l,exclude_neg=True,log_revd=True)

        afs.load_all_candidates()
        afs.count()
        q = afs.default_affine_decomp()
        sol = [((2, 0), {1, 2, 3, 5, 6}), ((3, 0), {8, 9, 10}), ((3, -33), {7}), ((3, -45), {4})]
        assert q == sol 

    def test__ModuloDecomp__afs_on_subsequence_(self):

        l = [2,5,11,13,14,29]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)

        aos = md.afs_on_subsequence_(0)
        sol1 = ((0, 6), [((2, 1), {1, 2, 5}), ((3, -25), {4}), ((3, -20), {3})])
        assert aos == sol1 

        # test for span format of `aos`
        q = aos[1]
        q2 = AffineFitSearch.decomp_to_span_fmt(q)
        sol2 = [[(2, 1), [1, 2]], [(3, -20), [3, 3]], [(3, -25), [4, 4]], [(2, 1), [5, 5]]]
        assert sol2 == q2 

    def test__ModuloDecomp__afs_on_subsequence(self):
        l = [2,5,11,4,14,44,6,27,3,15]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)

        md.afs_on_subsequence(0)
        md.afs_on_subsequence(1)
        md.afs_on_subsequence(2)
        md.afs_on_subsequence(3)

        sol = [((0, 3), [[(2, 1), [1, 2]]]),\
            ((3, 6), [[(3, 2), [4, 5]]]),\
            ((6, 8), [[(5, -3), [7, 7]]]),\
            ((8, 10), [[(5, 0), [9, 9]]])]

        assert md.afs_prt == sol 
        assert md.afs_prt_mod == [np.int32(19), np.int32(128), np.int32(129)]
        assert md.gleqvec_prt == [2, 5, 7, 9]

    def test__ModuloDecomp__continuous_merge(self): 
        
        # case 1 
        l = [6,0,-6]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)

        q = md.continuous_merge(False)
        assert md.gleqvec_prt == [2]
        assert md.afs_prt == [((0, 3), [[(1, 0), [1, 1]], [(2, -6), [2, 2]]])]
        assert md.afs_prt_mod == []

        # case 2 
        l = [6,0]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        q = md.continuous_merge(False)
        assert md.gleqvec_prt == [1]
        assert md.afs_prt == [((0, 2), [[(1, 0), [1, 1]]])]
        assert md.afs_prt_mod == []

        # case 3 
        l = [0,6]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        q = md.continuous_merge(False)
        assert md.gleqvec_prt == [1]
        assert md.afs_prt == [((0, 2), [[(1, 6), [1, 1]]])]
        assert md.afs_prt_mod == [] 

        # case 4 
        l = [6,0,6]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        q = md.continuous_merge(False)
        assert md.gleqvec_prt == [1, 2]
        assert md.afs_prt == [((0, 2), [[(1, 0), [1, 1]]]), ((2, 3), [])]
        assert md.afs_prt_mod == [np.int32(-6)]

        # case 5 
        l = [-6,0,6]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        q = md.continuous_merge(False)
        assert md.gleqvec_prt == [2]
        assert md.afs_prt == [((0, 3), [[(1, 0), [1, 1]], [(2, 6), [2, 2]]])]
        assert md.afs_prt_mod == []

    def test__ModuloDecomp__premerge_contiguous(self):
        # case: negative multiple
        l = [3,-7,23,-67]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)

        q = md.continuous_merge(False)
        assert md.gleqvec_prt == [1, 3]
        assert md.afs_prt_mod == [np.int32(-50)]
        assert md.afs_prt == [((0, 2), [[(2, -13), [1, 1]]]), ((2, 4), [[(2, -113), [3, 3]]])]

        md.premerge_contiguous(1,False)
        assert md.gleqvec_prt == [3]
        assert md.afs_prt_mod == []
        assert md.afs_prt == [((0, 4), [[(-3, 2), [1, 3]]])]

        # case: negative multiple #2 
        l = [3,-7,23,-67,203,-607,1]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)
        assert md.gleqvec_prt == [5, 6]
        assert md.afs_prt_mod == [np.int32(1822)]
        assert md.afs_prt == [((0, 6), [[(-3, 2), [1, 5]]]), ((6, 7), [])]

    def test__ModuloDecompRepr__reconstruct(self):
        l = [2,5,11,4,14,44,6,27,3,15]
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)

        mdr = ModuloDecompRepr(md)
        r = mdr.reconstruct()
        assert l == r 
        return

if __name__ == '__main__':
    unittest.main()
