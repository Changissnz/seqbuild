from desi.nvec_gen import * 
from intigers.extraneous import * 
from morebs2.numerical_generator import prg__constant,prg__n_ary_alternator 
import unittest


### lone file test 
"""
python3 -m tests.test_nvec_gen 
"""
###
class PointSetGen__TypeAffineMethods(unittest.TestCase):

    def test__PointSetGen__TypeAffine_generate_points__case1(self): 

        num_points = 12 
        ro_prg_ = prg__LCG(3,-13,42,990)
        ro_prg_ = prg__LCG(4,-13,42,990)
        ro_prg = prg__single_to_range_outputter(ro_prg_) 

        ro_prg2 = prg__LCG(37,-131,4211,9800)
        ro_prg2 = prg__single_to_range_outputter(ro_prg2) 


        prg = prg__constant(5000)
        prg = prg__LCG(6,3,2,350)
        psg = PointSetGen__TypeAffine(num_points,prg,ro_prg,ro_prg2)
        psg.generate_points(is_ordered=True,clear_data=True)

        sx = set() 
        for p in psg.point_seq: 
            sx |= {str(p)}
        assert len(sx) == len(psg.point_seq) 

        psg.set_num_points(9)
        psg.generate_points(is_ordered=True,clear_data=False) 


        prg = prg__n_ary_alternator(30,300,33)
        def prg_q():
            return prg() * prg() 
        psg.prg = prg_q 

        psg.set_new_gen_ad_parameters() 
        psg.set_new_gen_ad_parameters() 
        psg.set_new_gen_ad_parameters() 

        psg.set_num_points(17) 
        psg.generate_points(is_ordered=True,clear_data=False) 

        q = np.array(psg.point_seq[-17:]) 
        assert q.shape == (17,5)

        q2 = np.array(psg.point_seq[:-17])
        assert q2.shape == (len(psg.point_seq) - 17,9) 
        assert len(psg.point_seq) - 17 == 21 

        qs = set()
        for s in psg.ad_indices_seq:
            qs |= set(s) 
            assert len(s) > 0 
        assert list(qs) == [i for i in range(38)]

        psg.generate_points(is_ordered=False,clear_data=False) 
        assert len(psg.ad_indices_seq[-1]) == 0 
        assert len(psg.ad_indices_seq[-3]) == 0   

    def test__PointSetGen__TypeAffine_generate_points__case2(self): 

        num_points = 12 
        ro_prg_ = prg__LCG(3,-13,42,990)
        ro_prg_ = prg__LCG(4,-13,42,990)
        ro_prg = prg__single_to_range_outputter(ro_prg_) 

        ro_prg2 = prg__LCG(37,-131,4211,9800)
        ro_prg2 = prg__single_to_range_outputter(ro_prg2) 

        prg = prg__constant(5000)
        prg = prg__LCG(6,3,2,350)
        psg = PointSetGen__TypeAffine(num_points,prg,ro_prg,ro_prg2)
        psg.generate_points(is_ordered=True,clear_data=True)

        psg.gen_ad_parameters[2] = psg.gen_ad_parameters[0][0]
        psg.generate_points(is_ordered=True,clear_data=False)

        qx = np.array(psg.input_seq[-12:]) 
        assert qx.shape == (12,9)  
    

if __name__ == '__main__':
    unittest.main()