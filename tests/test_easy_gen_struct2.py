from .tbclf_gen_test_cases import * 
import unittest

### lone file test 
"""
py -m tests.test_easy_gen_struct2
"""
###
class TimeBasedCommLangFileGeneratorMethodsNumber2(unittest.TestCase):

    def test__TimeBasedCommLangFileGenerator__generate_CL_n2m__case1(self): 
        return 
        filepath = "test_auto.txt"
        use_prng_for_prng_pr = 0.5 

        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,use_prng_for_prng_pr)
        s,g = cauto.generate_CL_n2m()
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 1000
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)
        cauto.close()

    def test__TimeBasedCommLangFileGenerator__generate_CL_gg__case1(self):

        filepath = "test_auto.txt"
        use_prng_for_prng_pr = 0.5 

        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,use_prng_for_prng_pr)
        s,g = cauto.generate_CL_gg()
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 1000
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)

    def test__TimeBasedCommLangFileGenerator__generate_CL_afs__case1(self):
        return 
        filepath = "test_auto.txt"
        use_prng_for_prng_pr = 0.5 

        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,use_prng_for_prng_pr)
        s,g = cauto.generate_CL_afs()
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 1000
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)

    def test__TimeBasedCommLangFileGenerator__generate_CL_fit22__case1(self):

        filepath = "test_auto.txt"
        use_prng_for_prng_pr = 0.5 

        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,use_prng_for_prng_pr)
        s,g = cauto.generate_CL_fit22()
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 1000
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)

    def test__TimeBasedCommLangFileGenerator__generate_CL_lps__case1(self):

        filepath = "test_auto.txt"
        use_prng_for_prng_pr = 0.5 

        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,use_prng_for_prng_pr)
        s,g = cauto.generate_CL_lps()
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 1000
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)


if __name__ == '__main__':
    unittest.main()