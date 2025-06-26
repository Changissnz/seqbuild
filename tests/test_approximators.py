from desi.approximators import * 
from intigers.mod_prng import prg__iterable
from morebs2.numerical_generator import prg__n_ary_alternator,prg__constant
import unittest

### lone file test 
"""
python -m tests.test_approximators 
"""
###
class ApproximatorsMethods(unittest.TestCase):

    def test__LPSValueOutputter__next__case1(self):

        lcg_px = prg__LCG(13,56,31,907)

        xr = prg__single_to_nvec(lcg_px,2)
        l_out = prg__LCG(43,21,82,505)

        range_outputter = prg__LCG(41.3131,50.6767,31.0121,5007.421541) 
        range_outputter = prg__single_to_range_outputter(range_outputter)

        adjustment_type = 1 

        lvo = LPSValueOutputter(xr,l_out,range_outputter,\
                adjustment_type)

        cx = []
        qx = set()  
        for i in range(500): 
            q = next(lvo)
            cx.append(q)
            qx |= {q}

        assert len(qx) == len(cx) 

    def test__LPSValueOutputter__next__case2(self):

        lcg_px = prg__LCG(13,56,31,907)

        xr = prg__single_to_nvec(lcg_px,2)
        l_out = prg__LCG(43,21,82,505)

        range_outputter = prg__single_to_range_outputter(l_out)

        adjustment_type = 1 

        lvo = LPSValueOutputter(xr,l_out,range_outputter,\
                adjustment_type)

        cx = []
        qx = set()  
        for i in range(500): 
            q = next(lvo)
            cx.append(q)
            qx |= {q}

        assert len(qx) == 326 

#-----------------------------------

    def test__Fit22ValueOutputter__next__case1(self): 

        lcg_p1 = prg__LCG(-12,13.1,51.0707,2210) 
        lcg_p2 = prg__LCG(32,12.75635,649,1400)

        def px():
            return (lcg_p1(),lcg_p2()) 

        l_out = prg__iterable([5,7,3])
        r_out = prg__constant((0.,5000.)) 
        b_out = prg__n_ary_alternator(0,2,0) 
        point_conn_type = 1 
        adjustment_type = 2 
        fvo = Fit22ValueOutputter(px,l_out,r_out,b_out,point_conn_type,adjustment_type)

        qr,qr2 = [],set()
        for _ in range(50): 
            x = next(fvo)
            qr.append(x) 
            qr2 |= {x}
        assert len(qr2) == 46 



if __name__ == '__main__':
    unittest.main()
