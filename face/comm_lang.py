"""
base file for interface language design and corresponding 
interpretive functions to process commands. 
"""
from seqgen.gg_gen import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","idforest","optri",\
             "rch"]

LANG_KEYTERMS = ["make","run","with","set","for","iter","write","to"]  

def MAKE_proc(splitstr_cmd): 
    assert splitstr_cmd[0] == "make"

    if splitstr_cmd[1] == "lcg":
        assert splitstr_cmd[2] == "with" 
        parameters = splitstr_cmd[3] 
        parameters = parameters.split(",")
        assert len(parameters) == 4 

        parameters = tuple([float(p) for p in parameters])
        return prg__LCG(*parameters) 
    return None 