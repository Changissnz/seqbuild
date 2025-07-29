from desi.multi_metric import * 
from morebs2.numerical_generator import prg__LCG,\
    prg__n_ary_alternator,prg__constant
import unittest

### lone file test 
"""
python3 -m tests.test_multi_metric 
"""
###
class MultiMetricMethods(unittest.TestCase):

    """
    case demonstrates improvement of replication 
    from <ModuloDecomp> to <ModuloDecompV2>. 
    """
    def test__MultiMetric__agv2_measures__ngrammer__case1(self):
        prg = prg__LCG(-13,1400,-52421,3245329670)
        L = [prg() for _ in range(10**5)] 
        ML = MultiMetric(L)

        r = ML.agv2_measures__ngrammer(0,10,\
                    5,set_frange=True) 
        rsol = np.array([\
            [1., 0.52866],\
            [0.87572, 0.43247],\
            [0.85781, 0.35785],\
            [0.45176, 0.23741],\
            [0.60033, 0.22798],\
            [0.85968, 0.39586],\
            [0.81562, 0.44854],\
            [0.95448, 0.55513],\
            [0.82432, 0.58931],\
            [0.97958, 0.58931]])

        assert np.all(r[:,[0,1]] == rsol)

    def test__MultiMetric__summarize__case1(self):

        prg = prg__LCG(-13,1400,-52421,3245329670)
        L = [prg() for _ in range(10**2)] 
        ML = MultiMetric(L[:100])
        q = ML.summarize(ngram=10)
        assert equal_iterables(q[0],np.array([0.8373328 , 0.3676976 , 0.40110277]))
        assert q[1] == (65, 70) 
        assert q[2] == np.float64(0.7)

        prg2 = prg__LCG(13,4,14,3245329670)
        L2 = [prg2() for _ in range(10**2)] 
        ML2 = MultiMetric(L2[:100])
        q2 = ML2.summarize(ngram=10)
        assert equal_iterables(q2[0],np.array([0.8043259 , 0.3383135 , 0.29560167])) 
        assert q2[1] == (49, 61)
        assert q2[2] == np.float64(0.7)

        prg3 = prg__n_ary_alternator(-500,500,-500) 
        L3 = [prg3() for _ in range(10**2)] 
        ML3 = MultiMetric(L3[:100])
        q3 = ML3.summarize(ngram=10)
        assert equal_iterables(q3[0],np.array([0.124363  , 0.0703733 , 0.03626311]))
        assert q3[1] == (1, 1) and q3[2] == 0.7

        prg4 = prg__constant(6) 
        L4 = [prg4() for _ in range(10**2)] 
        ML4 = MultiMetric(L4[:100])
        q4 = ML4.summarize(ngram=10)
        assert q4[1] == (1,1) 
        assert q4[2] == 0.0 

        prg5 = prg__LCG(0,1,2,10000) 
        L5 = [prg5() for _ in range(10**2)] 
        ML5 = MultiMetric(L5[:100])
        q5 = ML5.summarize(ngram=10)
        assert q5[1] == (1,1) 
        assert q5[2] == 0.7

    def test__MultiMetric__load_mc_map__case1(self):

        prg6 = prg__LCG(2,2,0,10000)
        L6 = [prg6() for _ in range(10**2)] 
        ML6 = MultiMetric(L6[:100])
        q6 = ML6.summarize(ngram=10)
        ML6.load_mc_map()

        assert ML6.modcomplex_map[2] == (1,1)
        assert ML6.modcomplex_map[1] == (1,1)
        assert ML6.modcomplex_map[4] == (2,2)
        assert ML6.modcomplex_map[8] == (3,3)

        prg7 = prg__LCG(4,11,6,10000)
        L7 = [prg7() for _ in range(10**2)] 
        ML7 = MultiMetric(L7[:100])
        q7 = ML7.summarize(ngram=10)
        ML7.load_mc_map()

        assert ML7.modcomplex_map[54] == (66,73)
        assert ML7.modcomplex_map[34] == (67,72)
        assert ML7.modcomplex_map[27] == (72,77)



if __name__ == '__main__':
    unittest.main()


