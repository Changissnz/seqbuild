from intigers.lcg_v3 import * 
import unittest

def LCGV3__sample_E(add_noise:bool): 
    prg = prg__LCG(131,53,30,9007)

    start = 131 
    m = 53 
    a = 30 
    n0 = 0
    n1 = 9007 

    super_range = [-5000,5000]

    P = LCGV3(start,m,a,n0,n1,prg,super_range,ternary_size_range=DEFAULT_LCGV3_TERNARY_SIZE_RANGE,\
        ternary_delta_timestamp_range = DEFAULT_LCGV3_TERNARY_DELTA_TIMESTAMP_RANGE,\
        add_noise=add_noise,verbose=False)
    return P 

def LCGV3__correct_sign_change_count(P): 

    s = P.s_  
    correct_sign_change = 0 
    X = [] 

    for _ in range(2000): 
        i = P.activated_ternary.index
        s2 = P.activated_ternary[i]
        r = next(P) 
        X.append(r) 
        q = to_trinary_relation(r,s) 
        s = r  
        correct_sign_change += (s2 == q) 

    return correct_sign_change,X


### lone file test 
"""
py -m tests.test_lcg_v3 
"""
###
class LCGV3Methods(unittest.TestCase):
    
    def test__LCGV3__next__case1(self): 
        P = LCGV3__sample_E(False) 
        correct_sign_change,X = LCGV3__correct_sign_change_count(P) 

        assert correct_sign_change == 1882 
        assert len(set(X)) == 1308, "got {}".format(len(set(X)))

    def test__LCGV3__next__case2(self): 
        P = LCGV3__sample_E(True) 
        correct_sign_change,X = LCGV3__correct_sign_change_count(P) 

        assert correct_sign_change == 750, "got {}".format(correct_sign_change) 
        assert len(set(X)) == 1741, "got {}".format(len(set(X)))

if __name__ == '__main__':
    unittest.main()