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
        ##gv = gleqvec(intsq.l)

        aos = md.afs_on_subsequence_(0)
        sol1 = ((0, 6), [((2, 1), {1, 2, 5}), ((3, -25), {4}), ((3, -20), {3})])
        assert aos == sol1 

        # test for span format of `aos`
        q = aos[1]
        q2 = AffineFitSearch.decomp_to_span_fmt(q)
        sol2 = [[(2, 1), [1, 2]], [(3, -20), [3, 3]], [(3, -25), [4, 4]], [(2, 1), [5, 5]]]
        assert sol2 == q2 


if __name__ == '__main__':
    unittest.main()
