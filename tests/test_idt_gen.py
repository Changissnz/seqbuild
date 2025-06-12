from .idf_sample_z import * 
import unittest

### lone file test 
"""
python -m tests.test_idt_gen 
"""
###
class IDecForestMethods(unittest.TestCase):

    """
    tests for crash-less run, length of <IntSeq> output.  

    NOTE: inconclusive test 
    """
    def test__IDecForest__one_new_IntSeq__case1(self):

        _,q = IDecForest_sampleZ()

        assert len(q[0]) == 60, "Q0: {}".format(len(q[0])) 
        assert len(q[1]) == 23, "Q1: {}".format(len(q[1]))  
        assert len(q[2]) == 11, "Q2: {}".format(len(q[2])) 
        return 


    """
    tests for crash-less run, length of integer output sequences 

    NOTE: inconclusive test 
    """
    def test__IDecForest__process_seq_at_tree__case1(self):
        idf,q = IDecForest_sampleZ()
        x3,x4,x5 = q 

        x1,x2 = idf.ST[0] 
        qx = idf.process_seq_at_tree__sequential(x2,x1)

        x1_ = IntSeq(x1[:len(x3)])
        qx2 = idf.process_seq_at_tree__iso_sequential(x2,x1_,x3)
        qx3 = idf.process_seq_at_tree__inflow(x2,x1)

        assert len(qx) == 90, "got {}".format(len(qx))
        assert len(qx2) == 207, "got {}".format(len(qx2))
        assert len(qx3) == 127, "got {}".format(len(qx3))

        prg3 = prg__LCG(31,83,311,1500)
        ql = list(x1.l)
        ql = prg_seqsort(ql,prg3)
        x1_ = IntSeq(ql)
        qx4 = idf.process_seq_at_tree__iso_sequential(x2,x1_,x1)
        assert len(qx4) == len(qx3)
        return 

    def test__IDecForest__process_seq_at_tree__case2(self): 
        idf,q = IDecForest_sampleZ()
        idf.one_tree() 
        x1,x2 = idf.ST[0] 
        qx1 = idf.process_seq_at_tree__splat(x2,x1)
        qx2 = idf.process_seq_at_tree__splat(x2,x1)
        qx3 = idf.process_seq_at_tree__splat(x2,x1)

        assert len(qx1) == 19 
        assert len(qx1) == len(qx2) 
        assert len(qx2) == len(qx3)
        return

    def test__IDecForest__next__case1(self): 

        idf,q = IDecForest_sampleZ()

        print("DISPLAY TEST")
        for _ in range(1000): 
            x = next(idf)
    

if __name__ == '__main__':
    unittest.main()