from mini_dm.affine_io_fit import * 
from morebs2.numerical_generator import * 
from intigers.extraneous import * 
from desi.nvec_gen import * 
import unittest

def sample_affine_dataset(): 
    num_points = 12 
    ro_prg_ = prg__LCG(3,-13,42,990)
    ro_prg_ = prg__LCG(4,-13,42,990)
    ro_prg = prg__single_to_range_outputter(ro_prg_) 

    ro_prg2 = prg__LCG(37,-131,4211,9800)
    ro_prg2 = prg__single_to_range_outputter(ro_prg2) 

    prg = prg__constant(5000)
    prg = prg__LCG(6,3,2,350)
    psg = PointSetGen__TypeAffine(num_points,prg,ro_prg,ro_prg2)

    psg.generate_points(is_ordered=True,clear_data=False)
    return np.array(psg.input_seq), np.array(psg.point_seq) 

### lone file test 
"""
py -m tests.test_affine_io_fit
"""
###
class IOAffineFitMethods(unittest.TestCase):


    def test__IOAffineFit__init_hypotheses__case_1(self): 

        I,O = sample_affine_dataset() 

        prg = prg__LCG(23,-113,442,1990)
        
        # case 1 
        iaf = IOAffineFit(I,O,prg)

        ma_dim = (0,9)  
        ma_order = 0 
        cv = (0.5,0.5) 

        iaf.init_hypotheses(ma_dim,ma_order,cv) 

        for i in range(12): 
            F = iaf.i2f_map[i] 

            q = F.fit(I[i]) 
            y = O[i]
            assert equal_iterables(q,y)

        # case 2 
        ma_dim2 = (9,9)  
        iaf2 = IOAffineFit(I,O,prg)
        iaf2.init_hypotheses(ma_dim2,ma_order,cv) 

        for i in range(12): 
            F = iaf2.i2f_map[i] 

            q = F.fit(I[i]) 
            y = O[i]
            assert equal_iterables(q,y)

        # case 3 
        ma_order = 1 
        iaf2.init_hypotheses(ma_dim2,ma_order,cv) 

        for i in range(12): 
            F = iaf2.i2f_map[i] 

            q = F.fit(I[i]) 
            y = O[i]
            assert equal_iterables(q,y)

    """
    identity map, 1-to-1 
    """
    def test__IOAffineFit__init_hypotheses__case_2(self): 
        I,_ = sample_affine_dataset() 

        prg = prg__LCG(23,-113,442,1990)
        
        # case 1 
        iaf = IOAffineFit(I,I,prg)

        ma_dim = (0,0)  
        ma_order = 0 
        cv = (0.5,0.5) 

        iaf.init_hypotheses(ma_dim,ma_order,cv) 

        for i in range(12): 
            F = iaf.i2f_map[i] 

            q = F.fit(I[i]) 
            y = I[i]

            assert q == y 

    """
    random contribution vectors, 1-to-1 
    """
    def test__IOAffineFit__init_hypotheses__case_3(self): 
        I,O = sample_affine_dataset() 

        prg = prg__LCG(23,-113,442,1990)
        
        # case 1 
        iaf = IOAffineFit(I,O,prg)

        ma_dim = (0,9)  
        ma_order = 0 
        cv = "random" 

        iaf.init_hypotheses(ma_dim,ma_order,cv) 

        for i in range(12): 
            F = iaf.i2f_map[i] 

            q = F.fit(I[i]) 
            y = O[i] 
            assert equal_iterables(q,y)

if __name__ == '__main__':
    unittest.main()