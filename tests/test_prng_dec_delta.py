from desi.prng_dec_delta import * 
import unittest

def run_PRNGDecimalDelta_on_float(f,\
    num_iter=50000,enable_base_delta=False): 

    pdd = PRNGDecimalDelta(f,50,enable_base_delta)   
    L = [] 
    for _ in range(num_iter): 
        L.append(next(pdd)) 

    d = defaultdict(int) 
    for l in L: 
        d[l] += 1 

    d = [(k,v) for k,v in d.items()] 
    d = sorted(d,key=lambda x: x[1]) 
    return L, d 

### lone file test 
"""
py -m tests.test_prng_dec_delta  
"""
###
class ProcessSeqMethods(unittest.TestCase):

    def test__PRNGDecimalDelta__next__case_1(self):

        f = 3.040641711229946

        # case 1: base float remains constant 
        L,d = run_PRNGDecimalDelta_on_float(f,enable_base_delta=False)  

        assert d[-5:] == [(4604, 650), (1, 2891), (2, 3144), \
            (1.9767441860465116, 4313), (0, 9200)]        
        assert len(set(L)) == 14286

        # case 2: enabling base delta yields more unique values per 50K 
        L,d = run_PRNGDecimalDelta_on_float(f,enable_base_delta=True)  

        assert d[-5:] == [(4, 276), (1, 828), (2, 2888), \
            (1.9767441860465116, 3682), (0, 14085)]
        
        assert len(set(L)) == 22046

    def test__PRNGDecimalDelta__next__case_2(self):

        f = 10.011011010101011011

        # case 1 
        L,d = run_PRNGDecimalDelta_on_float(f,enable_base_delta=True)  

        assert d[-5:] == [(4, 275), (1, 832), (2.0, 2533), (1.9767441860465116, 3922), (0.0, 17769)]
        assert len(set(L)) == 18220

        # case 2 
        L,d = run_PRNGDecimalDelta_on_float(f,enable_base_delta=False)  

        assert d[-5:] == [(-1, 1765), (1.9767441860465116, 4080), (2.0, 5031), (1, 9229), (0.0, 22460)]
        assert len(set(L)) == 1484

if __name__ == '__main__':
    unittest.main()