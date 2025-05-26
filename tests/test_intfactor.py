from intigers.intfactor import * 
import unittest

def check_factor_for_one(ms,f): 
    c = 0
    for x in ms: 
        if f in x: 
            c += 1
        if c > 1: return False 

    return c != 0 

### lone file test 
"""
python -m tests.test_intfactor
"""
###
class IntFactorMethods(unittest.TestCase):

    def test__intersection_disjunction__seq_factors__case1(self):
        iseq3 = IntSeq([12,120,1200,12000]) 

        ms3 = factors_of_seq(iseq3)
        idsm3 = intersection_disjunction__seq_factors(ms3)
        assert idsm3[0] == {1, 2, 3, 4, 6, 12}
        for x in idsm3[1]: 
            assert check_factor_for_one(ms3,x) 

    def test__intersection_disjunction__seq_factors__case2(self):
        iseq2 = IntSeq([1242,23326,3248,6734634]) 
        ms2 = factors_of_seq(iseq2)
        idsm2 = intersection_disjunction__seq_factors(ms2)
        assert idsm2[0] == {1, 2}
        for x in idsm2[1]: 
            assert check_factor_for_one(ms2,x) 

    def test__intersection_disjunction__seq_factors__case3(self):

        iseq = IntSeq([100,120,400,360,600]) 
        ms = factors_of_seq(iseq)
        idsm = intersection_disjunction__seq_factors(ms)

        assert idsm[0] == {1, 2, 4, 5, 10, 20}
        for x in idsm[1]: 
            assert check_factor_for_one(ms,x) 

    def test__ISFactorSetOps__factor_count__case1(self): 
        L = [100,120,400,360,600]
        iseq = IntSeq(L) 
        ms = factors_of_seq(iseq)
        idsm = intersection_disjunction__seq_factors(ms)

        isfso = ISFactorSetOps(L)
        isfso.factor_count_() 

        assert isfso.cofactor_degree_set(1) == \
            {36, 72, 9, 360, 75, 300, 45, 400, 16, 80, 18, 180, 150, 600, 90}
        assert isfso.cofactor_degree_set(4) == {40, 8}
        assert isfso.cofactor_degree_set(5) == {1, 2, 4, 5, 20, 10}

        assert isfso.cofactor_degree_set(1) == idsm[1] 
        assert isfso.cofactor_degree_set(5) == idsm[0] 

    def test__ISFactorSetOps__primes__case1(self): 
        L = [5,7,0,1200,1,12000]
        isfso = ISFactorSetOps(L)
        isfso.factor_count_() 
        assert isfso.primes() == {5,7,0,1}
        assert type(isfso.is_prime(213123123123121)) == type(None) 

if __name__ == '__main__':
    unittest.main()