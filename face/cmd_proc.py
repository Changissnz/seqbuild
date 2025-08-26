
"""
code in this file mainly consists of wrapper functions around 
the seqbuild program's main functional features. 
"""

from .make_cmd import * 
from morebs2.matrix_methods import cr 
from mini_dm.nsfr import * 

import io 
import pickle 

LANG_KEYTERMS = ["make","run","with","set","for","iter","write",\
        "merge","to","open","convert","show","help","file","seq",\
        "obj","qualtest"]   

GENFORM_CONVERT_TYPES = ["range","ndim","nvec","tvec"]

def parse_object_to_str(O): 
    s = str(type(O)) 

    i0 = s.index("'") 
    s = s[i0+1:] 

    i1 = s.index("'") 
    s = s[:i1] 


    stat = "." in s 
    while stat: 
        i0 = s.index(".") 
        s = s[i0+1:] 
        stat = "." in s
    return s 

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
open file <FILENAME.txt>
open file <FILENAME.txt> for seq 
open file <FILENAME.txt> for obj
"""
def OPEN_proc(splitstr_cmd):
    assert splitstr_cmd[0] == "open" 
    assert splitstr_cmd[1] == "file" 

    # check for directory name
    dirPath = os.path.dirname(splitstr_cmd[2])

    write_modes = ["wb","wb"] # "wa"

    if len(splitstr_cmd) != 3: 
        assert len(splitstr_cmd) == 5
        assert splitstr_cmd[3] == "for" 
        assert splitstr_cmd[4] in {"seq","obj"} 
        
        if splitstr_cmd[4] == "seq": 
            write_modes = ["w","w"] # "a" 

    if not os.path.isdir(dirPath) and dirPath != "": 
        os.mkdir(dirPath)

    if not os.path.exists(splitstr_cmd[2]): 
        q = open(splitstr_cmd[2],write_modes[0])
        q.close() 

    return open(splitstr_cmd[2],write_modes[1]) 

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


def WRITE_seq(fi_obj,vector_length,seq,rounding_depth=0):
    s = ""
    
    if rounding_depth == 0: 
        castFunc = int 
    else: 
        assert rounding_depth > 0 and type(rounding_depth) == int 
        castFunc = round(float(x),rounding_depth)

    while len(seq) > 0:
        q = seq[:vector_length] 
        seq = seq[vector_length:] 
        s = vector_to_string(q,castFunc)
        fi_obj.write(s + "\n") 
    fi_obj.flush() 

"""
write <object> to <file_object>.
write <object> for <positive integer> iter to <file_object>. 
"""
def WRITE_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "write" 
    assert splitstr_cmd[1] in var_map
    var_obj = var_map[splitstr_cmd[1]] 

    assert splitstr_cmd[2] in {"to","for"}  
    
    if splitstr_cmd[2] == "to": 
        assert splitstr_cmd[3] in var_map 
        fi_obj = var_map[splitstr_cmd[3]]

        is_seq = type(var_obj) == list 
        assert isinstance(fi_obj,io.TextIOWrapper) if is_seq \
            else isinstance(fi_obj,io.BufferedWriter) 
    
        if is_seq: 
            WRITE_seq(fi_obj,10,var_obj,rounding_depth=0)
        else: 
            pickle.dump(var_obj,fi_obj)
    else: 
        num_iter = int(splitstr_cmd[3]) 
        assert num_iter > 0 

        assert splitstr_cmd[4] == "iter" 
        assert splitstr_cmd[5] == "to" 
        assert splitstr_cmd[6] in var_map 

        var_obj = MAIN_method_for_object(var_obj) 
        fi_obj = var_map[splitstr_cmd[6]]

        while num_iter > 0: 
            elements_per_line = min([10,num_iter]) 
            nx = [str(var_obj()) for _ in range(elements_per_line)]
            nx = ','.join(nx) 
            fi_obj.write(nx+"\n")
            num_iter -= elements_per_line

    fi_obj.flush() 
    fi_obj.close()    

def MERGE_proc(splitstr_cmd,var_map):

    assert splitstr_cmd[0] == "merge"
    
    generators = splitstr_cmd[1].split(",") 
    assert len(generators) >= 2 
    for (i,g) in enumerate(generators): 
        assert g in var_map 
        generators[i] = var_map[g] 
    
    assert splitstr_cmd[2] == "with" 

    operators = splitstr_cmd[3].split(",") 
    assert len(operators) == len(generators) - 1 
    ops = []
    for o in operators: 
        stat0 = o in ARITHMETIC_OP_STR_MAP
        stat1 = o in var_map 
        assert stat0 or stat1 

        if stat0: 
            ops.append(ARITHMETIC_OP_STR_MAP[o]) 
        else: 
            ops.append(var_map[o]) 
    
    def amalgamaschwann(): 
        i = 1 
        q = generators[0]() 
        while i < len(generators): 
            q2 = generators[i]() 
            op = ops[i-1]
            q = op(q,q2) 
            i += 1 
        return q 
    
    return amalgamaschwann

def LOAD_proc(splitstr_cmd): 
    assert splitstr_cmd[0] == "load" 
    assert len(splitstr_cmd) == 2 
    assert os.path.exists(splitstr_cmd[1]) 

    fi_obj = open(splitstr_cmd[1],'r') 
    rx = fi_obj.readlines()
    fi_obj.close() 
    return rx 

def QUALTEST_proc(splitstr_cmd,var_map):
    assert splitstr_cmd[0] == "qualtest" 

    if "," in splitstr_cmd[1]: 
        q = splitstr_cmd[1].split(",") 
        assert len(q) == 2
        assert q[0] in var_map 
        assert q[1] in var_map 

        G,G2 = var_map[q[0]],var_map[q[1]]
        G,G2 = MAIN_method_for_object(G),MAIN_method_for_object(G2)  
    else:
        assert splitstr_cmd[1] in var_map
        G = MAIN_method_for_object(var_map[splitstr_cmd[1]])
        G2 = None 

    assert type(G) in {MethodType,FunctionType}

    assert splitstr_cmd[2] == "for"  

    num_iter = int(splitstr_cmd[3]) 
    assert splitstr_cmd[4] == "iter" 

    gauge_depth = 1 
    deg_vec = None 
    if len(splitstr_cmd) > 5: 
        assert splitstr_cmd[5] == "with"
        assert len(splitstr_cmd) == 7 
        stat = "," in splitstr_cmd[6] 
        if stat: 
            deg_vec  = splitstr_cmd[6].split(",") 
            deg_vec = [int(d) for d in deg_vec] 
            gauge_depth = None
        else: 
            gauge_depth = int(splitstr_cmd[6])
            deg_vec = None 

    s = "COV,DIST,ENT,DIFF1,DIFF2,MCS\n"   
    m = 1  
    if type(G2) == type(None): 
        S = gauge_generator__MultiMetric(G,num_iter,gauge_depth,deg_vec,\
            set_frange=True,condense_ngram_output=True) 
        q = list(S[0][0])
        q.extend(S[0][1])
        q.append(S[0][2]) 
        q = np.array(q) 
        s += str(np.round(q)) + "\n"
    else: 
        S = cmp_generators__MultiMetric(G,G2,num_iter,gauge_depth,deg_vec,\
            set_frange=True)
        s += str(np.round(S[0],5)) + "\n"

        m = 2 

    # get the top 5 factors 
    if type(S[1]) != type(None): 
        q = [(k,np.mean(v)) for k,v  in S[1].items() if k not in {1,2}]  
        q = sorted(q,key=lambda x:x[1])[:5] 
        s += "\nmodular characteristic\n\n"
        for q_ in q: 
            s += str(q_[0]) + "\t" + str(np.array(S[1][q_[0]])) + "\n"
    
    # get the frequency map 
    s += "\nnumber of unique elements: " + str(len(S[2])) + "\n"
    s += "\ntotal number of elements: " + str(m * num_iter) + "\n" 
    return s

def CHAINTEST_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "chaintest" 

    G2 = None 
    # case: comparative test
    if "," in splitstr_cmd[1]: 
        q = splitstr_cmd[1].split(",") 
        assert len(q) == 2 
        assert q[0] in var_map and q[1] in var_map 
        G = MAIN_method_for_object(var_map[q[0]])
        G2 = MAIN_method_for_object(var_map[q[1]])
    # case: lone test 
    else: 
        assert splitstr_cmd[1] in var_map
        G = MAIN_method_for_object(var_map[splitstr_cmd[1]]) 

    assert splitstr_cmd[2] == "with"
    assert len(splitstr_cmd) == 4 

    qx = splitstr_cmd[3].split(',')
    assert len(qx) == 2 
    num_iter = int(qx[0])
    chain_length = int(qx[1]) 

    ag2 = APRNGGaugeV2(G,(0.,1.),0.5) 
    q2 = ag2.chaintest(num_iter,chain_length) 

    if type(G2) != type(None): 
        ag2 = APRNGGaugeV2(G2,(0.,1.),0.5)  
        q2_ = ag2.chaintest(num_iter,chain_length)
        q2 = q2 - q2_ 

    return q2 

"""
reads a file containing a numerical sequence 
--- 

read <filepath>. 
read <filepath> for <positive integer> iter.  
"""
def READ_proc(splitstr_cmd):
    assert splitstr_cmd[0] == "read" 
    assert os.path.exists(splitstr_cmd[1]) 

    l = len(splitstr_cmd) 
    assert l in {2,5}

    num_iter = 10 ** 5 
    if l == 5: 
        assert splitstr_cmd[2] == "for" 
        assert splitstr_cmd[4] == "iter" 
        num_iter = int(splitstr_cmd[3])
        assert num_iter > 0 

    fi_obj = open(splitstr_cmd[1],"r")
    nsfr = NSFileReader(fi_obj,float,False)

    qs = []
    for _ in range(num_iter): 
        q = next(nsfr) 
        if type(q) == type(None): break 

        qs.append(q) 
    nsfr.close() 
    return qs

"""
encrypt <filepath> with <outfile,prg>.
"""
def ENCRYPT_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "encrypt"
    assert os.path.exists(splitstr_cmd[1])
    assert splitstr_cmd[2] == "with" 

    qx = splitstr_cmd[3].split(",") 
    
    inf1 = splitstr_cmd[1] 
    outf1 = qx[0] 
    assert qx[1] in var_map 

    G = MAIN_method_for_object(var_map[qx[1]]) 
    sbc = SBCrypt(inf1,outf1,G,0.5) 
    while not sbc.fin_stat:
        sbc.encrypt_one_chunk()
    sbc.close() 