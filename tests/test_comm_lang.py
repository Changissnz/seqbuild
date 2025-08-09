from face.comm_lang import * 
import unittest

### lone file test 
"""
python -m tests.test_comm_lang 
"""
###
class CommLangMethods(unittest.TestCase):

    def test__MAKE_proc__case1(self):
        cmd = "make lcg with 400,532,31,4577"
        splitstr_cmd = cmd.split(" ")

        prg = MAKE_proc(splitstr_cmd)

        lx = []

        for _ in range(12): 
            lx.append(prg()) 

        assert equal_iterables(lx, \
            [400.0, 2289.0, 297.0, 2417.0, 4315.0, \
            2534.0, 2481.0, 1747.0, 304.0, 1564.0, \
            3642.0, 1504.0])

if __name__ == '__main__':
    unittest.main()
