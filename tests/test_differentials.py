from desi.differentials import * 
import unittest

### lone file test 
"""
python -m tests.test_differentials
"""
###
class DifferentialsMethods(unittest.TestCase):

    def test__PIDValueOutputter__next__case1(self):


        lcg_px = prg__LCG(131,53,30,9007)
        px = prg__single_to_nvec(lcg_px,2) 

        f_out = prg__LCG(400,12,17,450) 

        l_out = prg__LCG(43,21,82,505)

        lcg_px2 = prg__LCG(4.0321,-12.33333333,25.151515,800.55779)
        r_out = prg__single_to_range_outputter(lcg_px2) 

        #r_out = prg__constant((0.,5005.))

        adjustment_type = 1

        #-----------------------------------------------

        ps = [px() for _ in range(15)]
        ps = np.array(ps) 

        fvec = np.arange(0,14)
        fvec = fvec % 2 

        #pic = PointInterpolationContainer(ps,[1,1,0,0]) 

        pvo = PIDValueOutputter(px,f_out,\
                l_out,r_out,adjustment_type)
        #pvo.set_next_value() 

        xr = set() 
        for _ in range(1000): 
            xr |= {next(pvo)}
            #print(next(pvo))
        assert len(xr) == 1000 


if __name__ == '__main__':
    unittest.main()