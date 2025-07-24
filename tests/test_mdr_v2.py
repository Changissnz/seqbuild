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


if __name__ == '__main__':
    unittest.main()