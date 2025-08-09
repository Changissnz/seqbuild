
from seqgen.gg_gen import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","idforest","optri",\
             "rch"]

LANG_KEYTERMS = ["make","run","with","set","for","iter","write","to",\
            "open"]  

# TODO: incomplete 
def MAIN_object_method(q):

    if type(q) in {MethodType,FunctionType}:
        def f():
            return q() 
        return f 

    if type(q) in {LCGV2,LCGV3}:
        return q.__next__ 
    return -1 

# TODO: incomplete 
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

def RUN_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "run"
    assert splitstr_cmd[1] in var_map 
    
    # get the object from the variable name     
    var_obj = var_map[splitstr_cmd[1]] 

    f = MAIN_object_method(var_obj) 

    if len(splitstr_cmd) < 3:
        return f()

    lx = []

    assert splitstr_cmd[2] == "for" 

    iterations = int(splitstr_cmd[3]) 
    assert iterations > 0 
    assert splitstr_cmd[4] == "iter" 

    for _ in range(iterations):
        lx.append(f())
    return lx 

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
