from face.easy_gen_struct import * 

def TimeBasedCommLangFileGenerator__case_BLACKOUT(filepath,use_prng_for_prng_pr:float):
    fpath = "face/sample_script/"
    vector_files = [fpath + "dummy_obj.txt"] 
    comm_lang_files = [fpath + "commond_six.txt"]

    base_gen_name = "rx"

    return TimeBasedCommLangFileGenerator(filepath,base_gen_name,\
        use_prng_for_prng_pr,vector_files,comm_lang_files)

def TimeBasedCommLangFileGenerator__case_NULLSRC(filepath,use_prng_for_prng_pr:float):
    fpath = "face/sample_script/"
    vector_files = [] 
    comm_lang_files = []

    base_gen_name = "rx"

    return TimeBasedCommLangFileGenerator(filepath,base_gen_name,\
        use_prng_for_prng_pr,vector_files,comm_lang_files)

def assertion_on_TimeBasedCommLangFileGenerator_output(cauto,l,i): 
    cauto.write_out_to_CL_file()

    v = "set TESTVEC = run {} for {} iter.".format(l,i) 
    cauto.update_CL_file([v]) 
    cauto.write_out_to_CL_file() 

    Q = cauto.clp.vartable["TESTVEC"]
    assert type(Q) == list and len(Q) == i, "GOT {}".format(Q)

    for q in Q: 
        assert is_number(q), "got {}".format(q)
        
    cauto.close() 
    print("****************** FILEPATH: ",cauto.filepath)
    clp = CommLangParser(cauto.filepath)
    clp.process_file() 

    assert not clp.cmd_errors, "got {}".format(clp.cmd_errors)

    clp.close() 
    return 

