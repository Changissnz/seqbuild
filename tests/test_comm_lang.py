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
        var_map = dict() 

        prg = MAKE_proc(splitstr_cmd,var_map)

        lx = []

        for _ in range(12): 
            lx.append(prg()) 

        assert equal_iterables(lx, \
            [400.0, 2289.0, 297.0, 2417.0, 4315.0, \
            2534.0, 2481.0, 1747.0, 304.0, 1564.0, \
            3642.0, 1504.0])

    def test__OPEN_proc__case1(self):

        fp = "face/sample_script/commondeth.txt"
        cmd2 = "open file {}".format(fp)
        splitstr_cmd2 = cmd2.split(" ")
        fi = OPEN_proc(splitstr_cmd2)
        assert type(fi) == io.BufferedWriter 
        fi.close()

    def test__CommLangParser__process_command__case1(self): 
        clp = CommLangParser("face/sample_script/commond_one.txt") 

        clp.load_next_command()
        #'set G = make lcg with 400,532,31,4577.'
        q = clp.process_command()
        assert q[0] == "G" 
        assert type(q[1]) in {MethodType,FunctionType}

        clp.close() 

if __name__ == '__main__':
    unittest.main()
