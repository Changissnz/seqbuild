from seqgen.cs_gen import * 
from face.prng_stats import * 
import unittest

### lone file test 
"""
py -m tests.test_cs_gen 
"""
###
class CongruenceShifterGenMethods(unittest.TestCase):

    """
    statistical tests comparing different PRNG inputs into 
    a <CongruenceShifterGen>. 
    """
    def test__CongruenceShifterGen__next__case1(self):

        prg = prg__LCG(5,1121,17,5342)
        super_range = [-1000,1000]
        csg = CongruenceShifterGen(deepcopy(prg),deepcopy(super_range))

        # [0] 18 
        # [1] 11 
        # [2] 890
        prg_seq = [prg__LCG(3,14,-71,190), prg__LCG(5,5,1,1000), prg__LCG(5,1121,17,5342)] 

        pios = PRNGIOStats(prg_seq,"prg",csg,[4,11],5)
        pios.run() 

        assert pios.input_cmp(0,1,{1,2,3,4,5,6,7}) == 1 
        assert pios.input_cmp(0,2,{1,2,3,4,5,6,7}) == -1 
        assert pios.input_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

        assert pios.output_cmp(0,1,{1,2,3,4,5,6,7}) == -1
        assert pios.output_cmp(0,2,{1,2,3,4,5,6,7}) == -1 
        assert pios.output_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

        assert pios.io_cmp(0,{1,2,3,4,5,6,7}) == -1 
        assert pios.io_cmp(1,{1,2,3,4,5,6,7}) == 0 
        assert pios.io_cmp(2,{1,2,3,4,5,6,7}) == -1 


if __name__ == '__main__':
    unittest.main()



##csg2 = CongruenceShifterGen(deepcopy(prg),super_range,allow_stretch_and_shrink=True)

'''
>>> pv = 281
>>> av = 218.1922999999
>>> modrange_for_congruence(pv,av,[212.4384599999994, 293.7192299999997])
[0, -62.8077000001]
>>>
'''

