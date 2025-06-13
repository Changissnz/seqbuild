from intigers.lcg_v2 import * 
import unittest

### lone file test 
"""
python -m tests.test_lcg_v2 
"""
###
class LCGV2Methods(unittest.TestCase):

    def test__LCGV2__basefunctions__case1(self):

        # case 1
        l2 = LCGV2(2,2,0,0,10,0) 

        d1 = l2.init_dir()
        d2 = l2.gen_dir() 
        assert d1 == 1 
        assert d2 == 1 

        cv = l2.coverage()
        assert cv == 0.4 

        cl = l2.cycle_length()
        assert cl == 4 

        # case 2 
        l2 = LCGV2(0,2,1,0,10,0) 

        d1 = l2.init_dir()
        d2 = l2.gen_dir() 
        assert d1 == 1 
        assert d2 == 1 

        cv = l2.coverage()
        assert cv == 0.45 

        cl = l2.cycle_length()
        assert cl == 5 

    def test__LCGV2__io_map_summary__case1(self):

        l2 = LCGV2(0,2,1,0,10,0) 
        l2.io_map()
        l2.io_map_partition()
        l2.io_map_summary()

        assert len(l2.cycle_descriptors) == 5
        
        assert l2.cycle_descriptors[0].is_closed() and \
        l2.cycle_descriptors[2].is_closed()
        
        assert not l2.cycle_descriptors[1].is_closed() and \
            not l2.cycle_descriptors[3].is_closed() and \
            not l2.cycle_descriptors[4].is_closed() 

        dheads = {} 
        dheads[0] = {0, 1, 3, 5, 7}
        dheads[1] = {5}
        dheads[2] = {9}
        dheads[3] = {3}
        dheads[4] = {7}

        for (i,ld) in enumerate(l2.cycle_descriptors):
            assert ld.d["sub-cycle"] == dheads[i]

    def test__LCGV2__io_map_summary__case2(self):
        l2 = LCGV2(0,1,2,0,10,0) 
        l2.io_map()
        l2.io_map_partition()
        l2.io_map_summary()

        assert l2.gd.components == \
        [[{0}, {2}, {4}, {6}, {8}], [{1}, {3}, {5}, {7}, {9}]]
        assert l2.cycle_descriptors[0].is_closed() 
        assert l2.cycle_descriptors[0].d["sub-cycle"] == {0, 2, 4, 6, 8}
        assert l2.cycle_descriptors[1].is_closed() 
        assert l2.cycle_descriptors[1].d["sub-cycle"] == {1, 3, 5, 7, 9}        

if __name__ == '__main__':
    unittest.main()