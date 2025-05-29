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

    def test__ISFactorSetOps__remove_seq_elements(self): 
        L = [3002,4370,10006,912,431]
        isfso = ISFactorSetOps(L,int_limit=DEFAULT_INT_MAX_THRESHOLD)
        isfso.factor_count_() 

        ei = isfso.iseq.element_indices([3002,431])
        dxq = deepcopy(isfso.factor_count)
        dx = isfso.factorcount_for_elementindices(ei)

        #print(str(isfso))
        x = [isfso.factors[ei_] for ei_ in ei]
        x = flatten_setseq(x) 
        isfso.remove_seq_elements(set([3002,431])) 
        #print("----")
        #print(str(isfso))

        xr1 = numberdict_op(dxq,dx,f=sub)
        xr2 = isfso.factor_count
        assert equal_intdicts(xr1,xr2)
        assert not equal_intdicts(dxq,xr2)

        assert set(isfso.iseq.l) == {4370,10006,912}
        assert len(isfso.iseq) == 3 
        assert len(isfso.factors) == 3 

    def test__ISFactorSetOps__coprimes_of_case1(self):

        L = [3002,4370,10006,912,431]
        isfso = ISFactorSetOps(L,int_limit=DEFAULT_INT_MAX_THRESHOLD)
        isfso.factor_count_() 

        ei = isfso.iseq.element_indices([3002,431])
        dxq = deepcopy(isfso.factor_count)
        dx = isfso.factorcount_for_elementindices(ei)

        ans = {3002: {np.int32(431)},\
            4370: {np.int32(431)},\
            10006: {np.int32(431)},\
            912: {np.int32(431)},\
            431: {np.int32(912), np.int32(3002), \
                np.int32(4370), np.int32(10006)}}

        for l in L: 
            assert ans[l] == isfso.coprimes_of(l)

    def test__ISFactorSetOps__coprimes_of_case2(self):

        L = [480,320,6400,1280,804,7,19,11]
        isfso = ISFactorSetOps(L,int_limit=DEFAULT_INT_MAX_THRESHOLD)
        isfso.factor_count_() 

        ei = isfso.iseq.element_indices([3002,431])
        dxq = deepcopy(isfso.factor_count)
        dx = isfso.factorcount_for_elementindices(ei)

        ans = {480: {np.int32(11), np.int32(19), np.int32(7)},\
            320: {np.int32(11), np.int32(19), np.int32(7)},\
            6400: {np.int32(11), np.int32(19), np.int32(7)},\
            1280: {np.int32(11), np.int32(19), np.int32(7)},\
            804: {np.int32(11), np.int32(19), np.int32(7)},\
            7: {np.int32(480), np.int32(6400), np.int32(320), \
                np.int32(1280), np.int32(804), np.int32(19), np.int32(11)},\
            19: {np.int32(480), np.int32(6400), np.int32(320), \
                np.int32(1280), np.int32(804), np.int32(7), np.int32(11)},\
            11: {np.int32(480), np.int32(6400), np.int32(320), \
                np.int32(1280), np.int32(804), np.int32(19), np.int32(7)}}

        for l in L: 
            assert ans[l] == isfso.coprimes_of(l)

    def test__ISFactorSetOps__dsort(self):
        L = [3002,4370,10006,912,431]
        L = [480,320,6400,1280,804,7,19,11]
        isfso = ISFactorSetOps(L,int_limit=DEFAULT_INT_MAX_THRESHOLD)
        isfso.factor_count_() 

        dx = isfso.dsort(pkeys = [320,480]) 
        assert dx == [[3, 1], [6, 1], [12, 1], [15, 1], \
            [24, 1], [30, 1], [48, 1], [60, 1], [320, 1], \
            [64, 1], [480, 1], [96, 1], [240, 1], [120, 1], \
            [1, 2], [2, 2], [4, 2], [5, 2], [8, 2], [10, 2], \
            [16, 2], [20, 2], [160, 2], [32, 2], [40, 2], [80, 2]]

        dx2 = isfso.dsort(pkeys=None) 
        assert dx2 == [[15, 1], [24, 1], [30, 1], [480, 1], \
            [96, 1], [240, 1], [48, 1], [120, 1], [60, 1], \
            [6400, 1], [3200, 1], [1600, 1], [200, 1], \
            [400, 1], [25, 1], [800, 1], [100, 1], [50, 1], \
            [804, 1], [134, 1], [67, 1], [201, 1], [268, 1], \
            [402, 1], [7, 1], [19, 1], [11, 1], [3, 2], [6, 2], \
            [12, 2], [1280, 2], [256, 2], [128, 2], [640, 2], \
            [320, 3], [64, 3], [5, 4], [8, 4], [10, 4], [80, 4], \
            [16, 4], [20, 4], [160, 4], [32, 4], [40, 4], [2, 5], \
            [4, 5], [1, 8]]

        dx3 = isfso.median_sort(pkeys=None,r=0.5,fullpair_sequence=False) 
        assert dx3 == [402, 7, 268, 19, 201, 11, 67, 3, \
            134, 6, 804, 12, 50, 1280, 100, 256, 800, \
            128, 25, 640, 400, 320, 200, 64, 1600, 5, \
            3200, 8, 6400, 10, 60, 80, 120, 16, 48, 20, \
            240, 160, 96, 32, 480, 40, 30, 2, 24, 4, 15, 1],"got {}".format(dx3)


if __name__ == '__main__':
    unittest.main()