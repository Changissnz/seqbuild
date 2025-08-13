
"""
code in this file mainly consists of wrapper functions around 
the seqbuild program's main functional features. 
"""

from .make_cmd import * 
import io 

LANG_KEYTERMS = ["make","run","with","set","for","iter","write","to",\
            "open","convert"]   

GENFORM_CONVERT_TYPES = ["range","ndim","nvec","tvec"]

"""
run object 
run object for INTEGER iter 
run object with PARAMETERS 
"""
def RUN_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "run"
    assert splitstr_cmd[1] in var_map 
    
    # get the object from the variable name     
    var_obj = var_map[splitstr_cmd[1]] 

    f = MAIN_method_for_object(var_obj) 

    if len(splitstr_cmd) < 3:
        return f()

    lx = []

    # case: keyword `with`
    if splitstr_cmd[2] == "with":
        return f(splitstr_cmd[3]) 

    # case: keyword `for`,`iter`. 
    assert splitstr_cmd[2] == "for" 

    iterations = int(splitstr_cmd[3]) 
    assert iterations > 0 
    assert splitstr_cmd[4] == "iter" 

    for _ in range(iterations):
        lx.append(f())
    return lx 

"""
open file FILENAME.txt 
open new file FILENAME.txt 
"""
def OPEN_proc(splitstr_cmd):
    assert splitstr_cmd[0] == "open" 
    assert splitstr_cmd[1] == "file" 

    # check for directory name
    dirPath = os.path.dirname(splitstr_cmd[2])

    if not os.path.isdir(dirPath) and dirPath != "": 
        os.mkdir(dirPath)

    if not os.path.exists(splitstr_cmd[2]): 
        q = open(splitstr_cmd[2],"wb")
        q.close() 

    return open(splitstr_cmd[2],"ab")

"""
convert GENERATOR to type_x ?with arg1,...,argN?.

convert G to range
convert G to ndim with 5,6,9,10 
convert G to nvec with 6 
convert G to tvec with 7
""" 
def CONVERT_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "convert"

    assert splitstr_cmd[1] in var_map

    struct = var_map[splitstr_cmd[1]] 

    if not type(struct) in {MethodType,FunctionType}: 
        struct = MAIN_method_for_object(struct)
    
    assert splitstr_cmd[2] == "to"
    assert splitstr_cmd[3] in GENFORM_CONVERT_TYPES,"GOT \n{}".format(splitstr_cmd)

    if splitstr_cmd[3] == "range":
        return prg__single_to_range_outputter(struct) 

    assert splitstr_cmd[4] == "with" 

    if splitstr_cmd[3] == "ndim": 
        dim = string_to_vector(splitstr_cmd[5], castFunc = int)
        return prg__single_to_ndim_index_outputter(struct,dim)

    l = string_to_vector(splitstr_cmd[5],castFunc = int) 
    assert len(l) == 1
    l = l[0] 

    if splitstr_cmd[3] == "nvec": 
        return prg__single_to_nvec(struct,l) 

    def struct_():
        return int(round(struct())) 

    return prg__single_to_trinary_vector(struct_,l)

"""
write object to file_object 
"""
def WRITE_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "write" 
    assert splitstr_cmd[1] in var_map
    var_obj = var_map[splitstr_cmd[1]] 
    assert splitstr_cmd[2] == "to" 
    
    assert splitstr_cmd[3] in var_map 
    fi_obj = var_map[splitstr_cmd[3]]
    assert isinstance(fi_obj,io.BufferedWriter) 
    #isinstance(fi_obj,io.TextIOWrapper)
    
    pickle.dump(var_obj,fi_obj)
    fi_obj.close()