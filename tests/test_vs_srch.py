from mini_dm.vs_srch import * 
from desi.nvec_gen import * 
from morebs2.numerical_generator import prg__LCG
import unittest

def pointset_sample_z():
    num_points = 96
    ro_prg = prg__constant((-3,1200))
    ro_prg2 = prg__constant((-31,1101))

    prg = prg__LCG(16,7,5,3220) 
    psg = PointSetGen__TypeAffine(num_points,prg,ro_prg,ro_prg2)
    psg.set_mod_output(False)

    for _ in range(3):
        psg.set_new_gen_ad_parameters() 

    psg.generate_points(is_ordered=True,clear_data=False) 
    return deepcopy(psg )


### lone file test 
"""
python -m tests.test_vs_srch 
"""
###
class VSSearchMethods(unittest.TestCase):

    '''
    tests for crash-less execution
    '''
    def test__VSSearch__preproc(self):
        psg = pointset_sample_z() 
        prg = prg__LCG(16,7,5,3220) 
        vs = VSSearch(psg.input_seq,psg.point_seq,\
            None,None,None,prg)         
        vs2 = deepcopy(vs) 

        vs.preproc() 
        ma_dim = [0,8] 
        ma_order = 0

        vs2.preproc_v2(ma_dim,ma_order)
        
        assert len(vs.mahm.mhm.info) == 77 
        assert len(vs2.mahm.mhm.info) == 96 

    def test__VSSearch__move_one_hyp__uc__case1(self):
        psg = pointset_sample_z() 

        prg = prg__LCG(16,7,5,3220) 
        vs = VSSearch(psg.input_seq,psg.point_seq,\
            None,None,None,prg,sol_maxsize=10) 

        ma_dim = [0,8] 
        ma_order = 0
        vs.preproc_v2(ma_dim,ma_order) 

        vs.initial_hypotheses() 
        r = vs.move_one_hyp__uc(unit=100.,err_type=2)

        assert len(vs.search_queue) == 193,"got {}".format(len(vs.search_queue))
        assert len(vs.soln) == 10 
        assert len(vs.n2mac.ftable) == 1000 

        q = set([v[1] for v in vs.soln])
        assert len(q) == len(vs.soln) - 2,"len {}".format(len(q))

    def test__VSSearch__move_one_hyp__prg_guided(self):
        
        psg = pointset_sample_z() 
        prg = prg__LCG(16,7,5,3220) 
        vs = VSSearch(psg.input_seq,psg.point_seq,\
            None,None,None,prg,sol_maxsize=10,is_bfs_queue=False) 
        vs.preproc() 
    
        vs.initial_hypotheses() 

        for i in range(5): 
            print("move {}".format(i))
            r = vs.move_one_hyp__prg_guided(unit=100.,err_type=2)
        assert len(vs.n2mac.ftable) == 96 
    

if __name__ == '__main__':
    unittest.main()

