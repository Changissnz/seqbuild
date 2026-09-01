from mini_dm.ngram_ext import * 
from morebs2.numerical_generator import prg__LCG
import unittest

### lone file test 
"""
py -m tests.test_ngram_ext 
"""
###
class SequencePrePostMapMethods(unittest.TestCase):

    def test__SequencePrePostMapMethods__key_to_seq__case1(self): 
        prg = prg__LCG(12,3,2,73) ## 54  
        S = [prg() for _ in range(100)] 

        pre_size = 3
        post_size = 4 

        spm = SequencePrePostMap(S,pre_size,post_size) 
        spm.process() 

        """
        spm.M['38.0,43.0,58.0']
        """ 

        f0 = spm.key_to_seq('38.0,43.0,58.0')
        f1 = spm.key_to_seq__valuefreq_form('38.0,43.0,58.0')
        f2 = spm.key_to_seq__valuefreq_form('38.0,43.0,58.0',label_suffixes=True)
        f3 = spm.key_to_seq__set_form('38.0,43.0,58.0')
        f4 = spm.key_to_seq__labelled('38.0,43.0,58.0')

        assert np.all(np.array(f0) == np.array([30, 19, 59, 33, 30, \
            19, 59, 33, 30, 19, 59, 33, 30, 19, 59, 33, 30,\
            19, 59, 33, 30, 19, 59, 33, 30, 19, 59, 33, 30,\
            19, 59, 33])) 

        assert f1 == ['30.0,19.0,59.0,33.0', 8], "got {}".format(f1)
        assert f2 == [0, 8]

        assert np.all(np.array(f3) == np.array([19, 30, 33, 59, \
            19, 30, 33, 59, 19, 30, 33, 59, 19, 30, 33, 59, 19,\
            30, 33, 59, 19, 30, 33, 59, 19, 30, 33, 59, 19, 30,\
            33, 59]))

        assert f4 == [0, 0, 0, 0, 0, 0, 0, 0]


if __name__ == '__main__':
    unittest.main()


