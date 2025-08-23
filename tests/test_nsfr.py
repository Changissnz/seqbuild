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
        fr = NSFileReader(f_obj,int,False) 

        rx = []
        while not fr.is_finished: 
            n = next(fr) 
            if fr.is_finished: continue 
            rx.append(n) 

        assert len(rx) == 20 
        assert min(rx) == 36 
        assert max(rx) == 958

        fr.close() 

    def test__NSFileReader__process_command__case2(self): 
        fs = "dummy_file2.txt" 
        f_obj = open(fs,'r')
        sng = NSFileReader(f_obj,float,is_periodic=True)
        q = [next(sng) for _ in range(30)]
        assert q[:10] == q[20:] 
        sng.close() 

        f_obj = open(fs,'r')
        sng2 = NSFileReader(f_obj,float,is_periodic=False) 
        q = [next(sng2) for _ in range(30)]
        assert q[20:] == [None] * 10 
        sng2.close()

if __name__ == '__main__':
    unittest.main()
