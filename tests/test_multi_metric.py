from desi.multi_metric import * 
from morebs2.numerical_generator import prg__LCG
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

        assert np.all(r == rsol)


if __name__ == '__main__':
    unittest.main()


