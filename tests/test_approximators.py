from desi.approximators import * 
from intigers.mod_prng import prg__iterable
from morebs2.numerical_generator import prg__n_ary_alternator,prg__constant
import unittest

def Fit22ValueOutputter__sample_1PRNG(prg,point_conn_type,adjustment_type): 

    prg = prg__single_to_nvec(prg,2) 

    length_outputter = None 
    range_outputter = None 
    bool_outputter = None 

    F = Fit22ValueOutputter(prg,length_outputter,range_outputter,bool_outputter,\
        point_conn_type,adjustment_type) 
    return F 

### lone file test 
"""
py -m tests.test_approximators 
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

        assert len(qx) == 338, "got {}".format(len(qx))

    def test__LPSValueOutputter__next__case3(self): 
        # case 1: integer LCG w/ output in [0,907]
        lcg_px = prg__LCG(13.55,56,31,907.44)

        lcg_px_ = prg__single_to_int(lcg_px)
        xr = prg__single_to_nvec(lcg_px_,2)
        l_out = None
        r_out = None 
        adjustment_type = 1 

        lvo = LPSValueOutputter(xr,l_out,r_out,\
                adjustment_type=1) 

        L = [] 
        for _ in range(1000): 
            L.append(next(lvo)) 
        assert len(set(L)) == 517, "got {}".format( len(set(L)))

        # case 2: float version of the PRNG in case 1 
        xr = prg__single_to_nvec(lcg_px,2)
        lvo = LPSValueOutputter(xr,l_out,r_out,\
                adjustment_type=1) 

        L = [] 
        for _ in range(1000): 
            L.append(next(lvo)) 

        assert len(set(L)) == 1000, "got {}".format( len(set(L)))

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
        assert len(qr2) == 46, "got {}".format(len(qr2))

    def test__Fit22ValueOutputter__next__case2(self): 
        
        # case 1 
        prg = prg__n_ary_alternator(-100,100,-100)
        F = Fit22ValueOutputter__sample_1PRNG(prg,point_conn_type= 1,adjustment_type= 1)
        L = [] 
        for _ in range(1000):
            L.append(next(F))
        assert len(set(L)) == 156  

        # case 2 
        F = Fit22ValueOutputter__sample_1PRNG(prg,point_conn_type= 1,adjustment_type=2)
        L = [] 
        for _ in range(1000):
            L.append(next(F))
        assert len(set(L)) == 661, "got {}".format(len(set(L))) 

        # case 3 
        F = Fit22ValueOutputter__sample_1PRNG(prg,point_conn_type= 2,adjustment_type=1) 
        L = [] 
        for _ in range(1000):
            L.append(next(F))
        assert len(set(L)) == 146, "got {}".format(len(set(L))) 

        # case 4 
        F = Fit22ValueOutputter__sample_1PRNG(prg,point_conn_type= 2,adjustment_type=2) 
        L = [] 
        for _ in range(1000):
            L.append(next(F))
        assert len(set(L)) == 740, "got {}".format(len(set(L))) 

    """
    constant case: outputs [25,0*]
    """
    def test__Fit22ValueOutputter__next__case3(self): 
        prg = prg__constant(13) 
        F = Fit22ValueOutputter__sample_1PRNG(prg,point_conn_type= 1,adjustment_type=1) 
        
        L = [] 
        for _ in range(1000):
            L.append(next(F))
        assert len(set(L)) == 2, "got {}".format(len(set(L))) 
        return

if __name__ == '__main__':
    unittest.main()