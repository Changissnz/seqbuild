from .tbclf_gen_test_cases import * 
import unittest

### lone file test 
"""
py -m tests.test_easy_gen_struct
"""
###
class TimeBasedCommLangFileGeneratorMethods(unittest.TestCase):

    def test__TimeBasedCommLangFileGenerator__generate_CL_LCGv3__case1(self):
        # CASE: valid 
        filepath = "test_auto.txt"
        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.) 
        s,g = cauto.generate_CL_LCGv3() 
        cauto.update_CL_file(s)
        cauto.write_out_to_CL_file() 

        assert not cauto.clp.cmd_errors, "got {}".format(s) 
        R = cauto.clp.object_list("sequence",True)
        R2 = [r[0] for r in R]
        assert R2 == ['rx_vec_0', 'V', 'MV1', 'MV2', 'MV11', 'MV12', 'MV21', 'MV22']
        assert R2 == cauto.available_vectors
        
        G = cauto.clp.object_list("generator",True)
        G2 = [r[0] for r in G]
        assert G2 == ['G', 'M', 'M2', 'rx_0', 'rx_1']

        cauto.close() 

        clp = CommLangParser(filepath)
        clp.process_cmdlines() 

        assert not clp.cmd_errors
        clp.close() 
        return 

    def test__TimeBasedCommLangFileGenerator__preprocess__case1(self):
        vector_files = [] 
        comm_lang_files = ["roary.txt"] 
        filepath = "test_auto2.txt"
        base_gen_name = "rx"
        use_prng_for_prng_pr = 0. 

        cauto = TimeBasedCommLangFileGenerator(filepath,base_gen_name,\
            use_prng_for_prng_pr,vector_files,comm_lang_files)

        assert cauto.erroneous_files == ['roary.txt']
        assert cauto.errors_exist()

        cauto.close()
        return 

    def test__TimeBasedCommLangFileGenerator__generate_CL_optri__case1(self): 
        
        # subcase 1 
        filepath = "test_auto.txt"
        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.0)
        s,g = cauto.generate_CL_optri() 

        cauto.update_CL_file(s) 
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,"rx_1",20)

        # subcase 2 
        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.5)
        s,g = cauto.generate_CL_optri() 

        cauto.update_CL_file(s) 
        l = cauto.clp.single_output_generator_list()[-1] 
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,25) 
        return 

    def test__TimeBasedCommLangFileGenerator__generate_CL_qval__case1(self):
        filepath = "test_auto.txt"
        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.5)
        s,g = cauto.generate_CL_qval() 

        cauto.update_CL_file(s) 
        l = cauto.clp.single_output_generator_list()[-1] 
        i = 220 
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)
        return 

    def test__TimeBasedCommLangFileGenerator__generate_CL_mdrgen__case1(self):

        filepath = "test_auto.txt"
        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.5)
        s,g = cauto.generate_CL_mdrgen() 
        cauto.update_CL_file(s) 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 735 
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)
        return 

    def test__TimeBasedCommLangFileGenerator__generate_CL_rch__case1(self): 

        filepath = "test_auto_X.txt"
        i = 31 

        
        cauto = None 
        for _ in range(10): 
            cauto = TimeBasedCommLangFileGenerator__case_NULLSRC(filepath,0.5)
            s,g = cauto.generate_CL_rch() 
            cauto.update_CL_file(s) 

            l = cauto.clp.single_output_generator_list()[-1] 
            assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)

    def test__TimeBasedCommLangFileGenerator__generate_CL_pid__case1(self):

        filepath = "test_auto_X.txt"

        L = [40] * 5 + [200] * 5 

        for l_ in L: 
            cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.5)
            s,g = cauto.generate_CL_pid() 
            cauto.update_CL_file(s) 
            
            l = cauto.clp.single_output_generator_list()[-1] 
            assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,l_) 

    def test__TimeBasedCommLangFileGenerator__generate_CL_ssino__case1(self):

        filepath = "test_auto_Y.txt"

        L = [4000,3000] 

        for l_ in L: 
            cauto = TimeBasedCommLangFileGenerator__case_NULLSRC(filepath,0.5)
            s,g = cauto.generate_CL_ssino() 
            cauto.update_CL_file(s) 
            
            l = cauto.clp.single_output_generator_list()[-1] 
            assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,l_) 

    def test__TimeBasedCommLangFileGenerator__generate_CL_idforest__case1(self): 

        filepath = "test_auto_Z.txt"

        L = [5000,5000] 

        for l_ in L: 
            cauto = TimeBasedCommLangFileGenerator__case_NULLSRC(filepath,0.5)
            print("TT: ",cauto.t) 
            s,g = cauto.generate_CL_idforest() 
            cauto.update_CL_file(s) 
            
            l = cauto.clp.single_output_generator_list()[-1] 
            assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,l_) 


    def test__TimeBasedCommLangFileGenerator__generate_CL_weighted_pwop__case1(self): 

        filepath = "test_auto_X.txt"
        cauto = TimeBasedCommLangFileGenerator__case_NULLSRC(filepath,0.5)

        s,g = cauto.generate_CL_ssino()
        s1,g1 = cauto.generate_CL_weighted_pwop() 

        cauto.update_CL_file(s) 
        cauto.update_CL_file(s1) 
        cauto.write_out_to_CL_file() 

        assert len(cauto.clp.non_generators) == 1

        cauto.close() 

    def test__TimeBasedCommLangFileGenerator__generate_CL_shadow__case1(self): 

        filepath = "test_auto_A.txt"
        cauto = TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,0.5)

        s,g = cauto.generate_CL_shadow() 
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        i = 1000 
        #L = cauto.clp.single_output_generator_list()[-1] 
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,g,i) 

    def test__TimeBasedCommLangFileGenerator__generate_CL_LCG_bunch__case1(self): 

        filepath = "test_auto.txt"
        cauto = TimeBasedCommLangFileGenerator__case_NULLSRC(filepath,0.5)
        s,g = cauto.generate_CL_LCG_bunch() 
        cauto.update_CL_file(s) 
        cauto.write_out_to_CL_file() 

        l = cauto.clp.single_output_generator_list()[-1] 
        i = 1000
        assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i)



if __name__ == '__main__':
    unittest.main()