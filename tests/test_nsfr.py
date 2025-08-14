from mini_dm.nsfr import * 
import unittest

### lone file test 
"""
python -m tests.test_nsfr 
"""
###
class NSFileReaderMethods(unittest.TestCase):

    def test__NSFileReader__process_command__case1(self): 
        f = "dummy_file.txt" 
        f_obj = open(f,"r")

        castFunc = lambda x: round(float(x),2) 
        fr = NSFileReader(f_obj,int) 

        rx = []
        while not fr.fin_stat: 
            n = next(fr) 
            if fr.fin_stat: continue 
            rx.append(n) 

        assert len(rx) == 20 
        assert min(rx) == 36 
        assert max(rx) == 958

        fr.close() 


if __name__ == '__main__':
    unittest.main()
