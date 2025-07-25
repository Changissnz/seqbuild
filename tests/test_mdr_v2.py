from intigers.mdr_v2 import * 
import unittest

### lone file test 
"""
python3 -m tests.test_mdr_v2 
"""
###
class ModuloDecompV2Methods(unittest.TestCase):

    """
    case demonstrates improvement of replication 
    from <ModuloDecomp> to <ModuloDecompV2>. 
    """
    def test__ModuloDecompV2__fix_subend__case1(self):
        
        l2 = [2, 59, 140, 113, 122, 119, 20, 53, 242, 179]

        intsq = IntSeq(l2) 
        md2 = ModuloDecompV2(intsq) 

        mdr2 = ModuloDecompRepr(md2,2) 
        r2 = mdr2.reconstruct()
        assert np.all(r2 == l2) 

        mdr2_ = ModuloDecompRepr(md2,1) 
        r2_ = mdr2_.reconstruct()
        assert not np.all(r2_ == l2) 
        
        return 

    def test__ModuloDecompV2__reconstruct__case1(self):

        # case: <ModuloDecomp> fits L 
        L = [12,-31,14,-2,10,100,112]

        md = ModuloDecomp(IntSeq(L))
        md.merge(False)
        mdr0 = ModuloDecompRepr(md,1)
        r0 = mdr0.reconstruct()
        assert np.all(r0 == L) 

        # case: <ModuloDecomp> does not fit L2 
        L2 = [-100,-50,-90,-80,-110,130] 
        md2 = ModuloDecomp(IntSeq(L2))
        md2.merge(False)

        mdr = ModuloDecompRepr(md2,1)
        r = mdr.reconstruct()
        assert not np.all(r == L2) 

        # case: <ModuloDecompV2> fits L3 and 
        #       <ModuloDecomp> does not. 
        L3 = [-79,-50,-90,-80,-110,130] 
        md3_ = ModuloDecompV2(IntSeq(L3))
        mdr3_ = ModuloDecompRepr(md3_,2)
        r3_ = mdr3_.reconstruct()
        assert np.all(r3_ == L3)

        md3 = ModuloDecomp(IntSeq(L3))
        md3.merge(False)
        mdr3 = ModuloDecompRepr(md3,1)
        r3 = mdr3.reconstruct()
        assert not np.all(r3 == L3) 
        return  


if __name__ == '__main__':
    unittest.main()