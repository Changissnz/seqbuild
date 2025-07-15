from mini_dm.puc import * 
import unittest

### lone file test 
"""
python -m tests.test_puc 
"""
###
class PUCrawlerMethods(unittest.TestCase):

    def test__PUCrawler__next__case1(self):

        v = np.zeros((12,)) 
        u = 0.05 
        partition = [2,4,6] 

        uc = PUCrawler(v,u,partition)

        q = next(uc) 

        sol1 = [np.array([-0.05, -0.05]),
                np.array([-0.05, -0.05, -0.05, -0.05]),\
                np.array([-0.05, -0.05, -0.05, -0.05, -0.05, -0.05])]

        for(i,q_) in enumerate(q):
            assert np.all(q_ == sol1[i]) 

        s = None 
        for i in range(2 ** 12):
            s_ = next(uc) 
            if type(s_) == type(None): 
                break 
            s = s_ 

        for(i,q_) in enumerate(s):
            assert np.all(q_ == -sol1[i]) 


if __name__ == '__main__':
    unittest.main()