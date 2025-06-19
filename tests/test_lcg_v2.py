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
        dheads[0] = {1, 3, 5, 7}
        dheads[1] = {2, 5}
        dheads[2] = {9}
        dheads[3] = {3, 6}
        dheads[4] = {8, 7}

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
        assert type(l2.cycle_descriptors[0].d["sub-cycle"]) == type(None)
        assert l2.cycle_descriptors[1].is_closed() 
        assert type(l2.cycle_descriptors[1].d["sub-cycle"]) == type(None)


    def test__APRNGGaugeV2__matrix_cat_entropy__case1(self): 

        mx2 = np.array([[160., 376., 120., 366.],\
            [ 80., 356.,  40., 346.],\
            [150., 336., 110., 326.],\
            [ 70., 316.,  30., 382.],\
            [140., 372., 100., 362.],\
            [ 60., 352.,  20., 342.],\
            [130., 332.,  90., 322.],\
            [ 50., 386., 160., 376.],\
            [120., 366.,  80., 356.],\
            [ 40., 346., 150., 336.],\
            [110., 326.,  70., 316.],\
            [ 30., 382., 140., 372.],\
            [100., 362.,  60., 352.],\
            [ 20., 342., 130., 332.]])

        franges = (0.,400.)
        is_rowwise = True 
        is_local_frange = True 
        sl_info = 5 
        count_type = "absdiff" 
        round_depth = 4 

        evec = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec == [0.6,0.8,0.6,0.6667,0.6,0.8,0.6,0.6667,\
            0.7333,0.6667,0.6667,0.6667,0.7333,0.6667])
        
        is_rowwise = False 
        evec2 = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec2 == np.array([0.5,0.0308, 0.4615, 0.0308]))

        franges = None 
        evec3 = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec3 == np.array([0.5077, 0.3846, 0.5231, 0.4308]))

        franges = (500.,1000.)
        evec4 = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec4 == np.array([0.2   , 0.    , 0.1846, 0.    ])) 

        franges = np.array([[0.,200.],[350.,400.],[0.,200.],[350.,400.]])
        evec5 = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec5 == np.array([0.4808, 0.3846, 0.4808, 0.3462])) 

        sl_info = np.array([1,3,6,6],dtype=np.int32)
        evec6 = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec6 == np.array([0.,0.2821, 0.4308, 0.3538]))

        franges = np.array([[0.,200.],[0.,400.],[0.,200.],[0.,400.]])
        evec7 = APRNGGaugeV2__matrix_cat_entropy(mx2,franges,is_rowwise,is_local_frange,\
            sl_info,count_type,round_depth)
        assert np.all(evec7 == np.array([0.    , 0.    , 0.4308, 0.0897]))


if __name__ == '__main__':
    unittest.main()