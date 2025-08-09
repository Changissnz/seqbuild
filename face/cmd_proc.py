
from seqgen.gg_gen import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","idforest","optri",\
             "rch"]

LANG_KEYTERMS = ["make","run","with","set","for","iter","write","to",\
            "open"]  

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

def RUN_proc(splitstr_cmd): 
    return -1  

def OPEN_proc(splitstr_cmd):
    assert splitstr_cmd[0] == "open" 
    assert splitstr_cmd[1] == "file" 

    # check for directory name
    dirPath = os.path.dirname(splitstr_cmd[2])
    if not os.path.isdir(dirPath): 
        os.mkdir(dirPath)

    if not os.path.exists(splitstr_cmd[2]): 
        open(splitstr_cmd[2],"w")

    return open(splitstr_cmd[2],"a")
