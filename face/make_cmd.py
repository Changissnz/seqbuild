from seqgen.gg_gen import * 
from seqgen.mdr_gen import * 
from seqgen.optri_gen import * 
from seqgen.rch_gen import * 
from seqgen.idt_gen import * 
from seqgen.shadow_gen import * 
from desi.fraction import * 
from desi.differentials import * 
from intigers.lcg_v3 import * 
from intigers.prng_pw_op import *
from seqgen.sb_crypt import * 
from seqgen.ssi_netop import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","mdrv2",\
        "mdrgen","iomaps","idforest","optri","rch","qval",\
        "pid","echo","shadow","ssino"]  

ARITHMETIC_OP_STR_MAP = {"+":add,"-":sub,"/":zero_div0,"*":mul,"%":mod} 

def is_stringized_number(n): 
    istat = False 
    try: 
        int(n)
        istat = True 
    except: 
        pass 
    return istat 

# TODO: incomplete 
def MAIN_method_for_object(q):

    if type(q) in {MethodType,FunctionType}:
        def f():
            return q() 
        return f 

    if type(q) in {LCGV2,LCGV3}:
        return q.__next__ 

    if type(q) == MultiMetric: 
        
        def f(ngram): 
            ngram_ = int(ngram) 
            assert ngram_ > 0 
            qx = q.summarize(ngram_,condense_ngram_output=True)
            
            try: q.load_mc_map()
            except: print("excessively large numbers --> cannot load MC map")

            return qx 
        return f

    if type(q) == ModuloDecompRepr:
        def f(first):
            q.reset_first(int(first),False)
            print("resetting")
            return q.reconstruct() 
        return f 

    if type(q) == MDRGen: 
        return q.__next__ 

    if type(q) == OpTriGen: 
        return q.__next__ 

    if type(q) == PIDValueOutputter: 
        return q.__next__ 

    if type(q) == QValueOutputter:
        return q.__next__ 

    if type(q) == RCHAccuGen:
        return q.__next__ 

    if type(q) == NSFileReader:
        return q.__next__ 

    if type(q) == ShadowGen:
        return q.__next__ 

    if type(q) == ModPRNGOutputter: 
        return q.__next__ 

    if type(q) == IDecForest: 
        return q.__next__

    if type(q) == SSINetOp: 
        return q.__next__

    return -1 

def MAKE_lcgvx(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make" 
    assert "lcg" in splitstr_cmd[1] 

    if splitstr_cmd[1] == "lcg":
        assert splitstr_cmd[2] == "with" 
        parameters = splitstr_cmd[3] 
        parameters = parameters.split(",")
        assert len(parameters) == 4 

        parameters = tuple([float(p) for p in parameters])
        return prg__LCG(*parameters) 

    if splitstr_cmd[1] == "lcgv2":
        assert splitstr_cmd[2] == "with" 

        parameters = splitstr_cmd[3] 
        parameters = parameters.split(",")
        
        assert len(parameters) == 5 or len(parameters) == 7

        if len(parameters) == 5:
            sc_size = 50 
            preproc_gd = False
        else: 
            sc_size = parameters.pop(5) 
            preproc_gd = parameters.pop(5) 

            sc_size = int(sc_size) 
            assert sc_size > 0 
            
            preproc_gd = bool(int(preproc_gd)) 
        
        parameters = [float(p) for p in parameters] 

        return LCGV2(parameters[0],parameters[1],parameters[2],\
            parameters[3],parameters[4],sc_size,preproc_gd=preproc_gd)

    if splitstr_cmd[1] == "lcgv3": 
        assert splitstr_cmd[2] == "with" 

        parameters = splitstr_cmd[3] 
        parameters = parameters.split(",")

        assert len(parameters) in {5,8,9,10},"parameters len {}: {}".format(len(parameters),parameters)
        px = [float(p) for p in parameters[:5]] 

        # case: basic case, functionally identically to standard LCG 
        if len(parameters) == 5: 
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],0,False,None,None,\
                False,False) 
            return g3 

        assert parameters[5] in var_map 
        prg = var_map[parameters[5]]
        
        try: 
            super_range = (float(parameters[6]),float(parameters[7])) 
            assert is_valid_range(super_range,False,False) 
        except: 
            super_range = None 

        # case: trinary guided
        if len(parameters) == 8: 
            exclude_zero__auto_td = True
            is_rmod = False 

            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],10,False,prg,super_range,\
                exclude_zero__auto_td,is_rmod)

        # case: trinary guided w/ option to exclude 0 from trinary vector guides. 
        elif len(parameters) == 9: 
            exclude_zero__auto_td = bool(int(parameters[8]))
            is_rmod = False 

            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],10,False,prg,super_range,\
                exclude_zero__auto_td,is_rmod)

        # case: trinary guided w/ option to exclude 0 from trinary vector guides, 
        # option to use `reflective modification` output mode 
        else: 
            exclude_zero__auto_td = bool(int(parameters[8]))
            is_rmod = bool(int(parameters[9]))
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],10,False,prg,super_range,\
                exclude_zero__auto_td,is_rmod)

        m0 = int(prg()) % 2 + 1 
        m1 = int(prg()) % 26 + 5 
        m2 = bool(int(prg()) % 2)
        g3.static_autoset(m0,m1,m2) 
        return g3     
        
    assert False 

def MAKE_mdrvx(splitstr_cmd,var_map): 
    if splitstr_cmd[1] == "mdr": 
        assert splitstr_cmd[2] == "with" 
        assert splitstr_cmd[3] in var_map 
        lx = var_map[splitstr_cmd[3]] 
        mdx = ModuloDecomp(IntSeq(lx)) 
        mdx.merge(False)
        return ModuloDecompRepr(mdx,reconstruct_type=1)
    
    if splitstr_cmd[1] == "mdrv2": 
        assert splitstr_cmd[2] == "with" 

        parameters = splitstr_cmd[3] 
        parameters = parameters.split(",")
        assert parameters[0] in var_map 

        exclude_neg = False 

        if len(parameters) == 2:
            i = bool(int(parameters[1]))
            exclude_neg = i 
        
        lx = var_map[parameters[0]]
        mdx = ModuloDecompV2(IntSeq(lx),exclude_neg)
        return ModuloDecompRepr(mdx,reconstruct_type=2) 

    if splitstr_cmd[1] == "mdrgen": 
        assert splitstr_cmd[2] == "with"

        parameters = splitstr_cmd[3] 
        parameters = parameters.split(",")
        assert parameters[0] in var_map 
        mdr = var_map[parameters[0]] 
        assert type(mdr) == ModuloDecompRepr 

        assert parameters[1] in var_map
        prg = var_map[parameters[1]]
        prg = MAIN_method_for_object(prg) 

        lp = len(parameters[2:])
        assert lp in {2,7} 
        exclude_neg = bool(int(parameters[2])) 
        gentype = int(parameters[3]) 
        assert gentype in {1,2} 

        gt2_rcswitch,gt2_sel1 = True,True
        gt2_sel2,gt2_sel3 = True,True
        gt2_seed_in_output = False  
        if lp == 7: 
            gt2_rcswitch = bool(int(parameters[4])) 
            gt2_sel1 = bool(int(parameters[5])) 
            gt2_sel2 = bool(int(parameters[6]))
            gt2_sel3 = bool(int(parameters[7])) 
            gt2_seed_in_output = bool(int(parameters[8]))

        return MDRGen(mdr,prg,exclude_neg,gentype,\
        gt2_rcswitch,gt2_sel1,gt2_sel2,\
        gt2_sel3,gt2_seed_in_output)

    assert False 

def MAKE_op2(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make"
    assert splitstr_cmd[1] == "op2" 
    assert splitstr_cmd[2] == "with" 

    parameters = splitstr_cmd[3].split(",") 
    assert len(parameters) in {1,3,4}

    if len(parameters) == 1:
        assert parameters[0] in ARITHMETIC_OP_STR_MAP,"NOO." 
        return ARITHMETIC_OP_STR_MAP[parameters[0]] 

    if parameters[0] in ARITHMETIC_OP_STR_MAP:
        op1 = ARITHMETIC_OP_STR_MAP[parameters[0]] 
    else: 
        assert parameters[0] in var_map,"wrongo"
        op1 = var_map[parameters[0]] 

    if parameters[1] in ARITHMETIC_OP_STR_MAP:
        op2 = ARITHMETIC_OP_STR_MAP[parameters[1]] 
    else: 
        assert parameters[1] in var_map,"wrongo #2"
        op2 = var_map[parameters[1]] 

    weight = float(parameters[2]) 
    order = 0 
    if len(parameters) == 4: 
        order = int(parameters[3]) 
        assert order in {0,1,2} 
    return one_weighted_pairwise_operator(op1,op2,weight,order)

"""
int,prg,float
num nodes,prg,mutation rate  

int,prg,float,int 
num nodes,prg,mutation rate,queue capacity 

int,prg,float,int,float,float 
num nodes,prg,mutation rate,queue capacity,output range 
"""
def MAKE_rch(splitstr_cmd,var_map):
    assert splitstr_cmd[0] == "make"
    assert splitstr_cmd[1] == "rch" 
    assert splitstr_cmd[2] == "with" 

    q = splitstr_cmd[3].split(",") 
    l = len(q) 
    assert l in {3,4,6} 

    num_nodes = int(q[0]) 
    
    assert q[1] in var_map 
    prg = var_map[q[1]] 
    prg = MAIN_method_for_object(prg) 
    mut_rate = float(q[2]) 

    queue_capacity = 1000 
    output_range = [-1000,1000]
    if l >= 4: 
        queue_capacity = int(q[3]) 
        assert queue_capacity > 2  
    
    if l == 6: 
        output_range[0] = float(q[4])
        output_range[1] = float(q[5])
        assert output_range[1] > output_range[0] 

    output_range = tuple(output_range)

    orange = (3,9)
    ufreq_range = (2,11) 
    rch = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,orange,prg,\
        ufreq_range,mut_rate,queue_capacity)
    rch.output_range = output_range
    return rch 

"""
G_1,G_2,...,G_j
""" 
def MAKE_iomaps(splitstr_cmd,var_map): 
    assert splitstr_cmd[1] == "iomaps" 
    assert splitstr_cmd[2] == "with" 

    X = splitstr_cmd[3].split(",") 
    assert len(X) > 0 

    G = [] 
    for x in X: 
        assert x in var_map 
        r = MAIN_method_for_object(var_map[x])
        rx = lambda x2: r() + x2
        G.append(rx)  

    mpo = ModPRNGOutputter(G)
    return mpo

"""
list,ModPRNGOutputter,int,int range,max number of trees,prg,prg2
list,ModPRNGOutputter,int,int range,max number of trees,prg

list,ModPRNGOutputter,int,int range,prg,prg2
list,ModPRNGOutputter,int,int range,prg

list,ModPRNGOutputter,int,max number of trees,prg,prg2
list,ModPRNGOutputter,int,max number of trees,prg

list,ModPRNGOutputter,int,prg,prg2
list,ModPRNGOutputter,int,prg

list,ModPRNGOutputter,prg,prg2
list,ModPRNGOutputter,prg
"""
def MAKE_idforest(splitstr_cmd,var_map): 
    assert splitstr_cmd[1] == "idforest" 
    assert splitstr_cmd[2] == "with" 

    parameters = splitstr_cmd[3].split(",") 
    L = len(parameters)

    assert L in {8,7,6,5,4,3}
    assert parameters[0] in var_map
    V = var_map[parameters[0]] 

    assert type(V) == list 
    assert parameters[1] in var_map
    MO = var_map[parameters[1]] 
    assert type(MO) == ModPRNGOutputter 

    idf = None 
    cache_size = 100 
    reprod_rate_range = (15,50)         
    max_trees = 10 
    G = None  
    prg2 = None 

    if L == 3: 
        assert parameters[2] in var_map
        G = var_map[parameters[2]] 
        assert type(G) in {MethodType,FunctionType}  
    elif L == 4:
        # case
        istat = is_stringized_number(parameters[2])

        if istat: 
            cache_size = int(parameters[2]) 
            assert cache_size >= 2

            G = var_map[parameters[3]] 
            assert type(G) in {MethodType,FunctionType}
        # case 
        else: 
            G = var_map[parameters[2]] 
            assert type(G) in {MethodType,FunctionType}  

            assert parameters[3] in var_map
            prg2 = var_map[parameters[3]] 
            assert type(prg2) in {MethodType,FunctionType}  

    elif L == 5:
        cache_size = int(parameters[2]) 

        istat = is_stringized_number(parameters[3]) 
        if istat: 
            max_trees = int(parameters[3])
            G = var_map[parameters[4]] 
            assert type(G) in {MethodType,FunctionType} 
        else: 
            G = var_map[parameters[3]] 
            assert type(G) in {MethodType,FunctionType} 

            prg2 = var_map[parameters[4]]
            assert type(prg2) in {MethodType,FunctionType} 
    elif L == 6:
        cache_size = int(parameters[2]) 
        x = int(parameters[3]) 

        istat = is_stringized_number(parameters[4]) 
        if istat: 
            x2 = int(parameters[4]) 
            reprod_rate_range = (x,x2) 
            assert reprod_rate_range[0] <= reprod_rate_range[1] 

            G = var_map[parameters[5]] 
            assert type(G) in {MethodType,FunctionType} 
        else: 
            max_trees = x 
            assert max_trees > 0 

            G = var_map[parameters[4]] 
            assert type(G) in {MethodType,FunctionType} 

            prg2 = var_map[parameters[5]] 
            assert type(G) in {MethodType,FunctionType} 
    elif L >= 7:
        cache_size = int(parameters[2]) 
        r0 = int(parameters[3]) 
        r1 = int(parameters[4]) 
        assert r0 < r1         
        reprod_rate_range = (r0,r1) 

        istat = is_stringized_number(parameters[5])
        if istat: 
            max_trees = int(parameters[5])
            G = var_map[parameters[6]] 
            assert type(G) in {MethodType,FunctionType} 

            if L == 8: 
                prg2 = var_map[parameters[6]] 
                assert type(prg2) in {MethodType,FunctionType} 
        else: 
            G = var_map[parameters[5]] 
            assert type(G) in {MethodType,FunctionType} 

            prg2 = var_map[parameters[6]] 
            assert type(prg2) in {MethodType,FunctionType} 

    idf = IDecForest(IntSeq(V),MO,cache_size,reprod_rate_range,max_trees,G,prg2,False)    
    return idf

def MAKE_echo(splitstr_cmd):
    assert splitstr_cmd[1] == "echo" 
    assert splitstr_cmd[2] == "with" 

    is_periodic = False 
    file_path = None 
    if "," in splitstr_cmd[3]:
        qx = splitstr_cmd[3].split(",") 
        assert len(qx) == 2 

        is_periodic = bool(int(qx[0]))
        file_path = qx[1] 
    else: 
        file_path = splitstr_cmd[3]

    assert os.path.exists(file_path) 
    fi_obj = open(file_path,'r')
    nsfr = NSFileReader(fi_obj,float,is_periodic)
    return nsfr 

def MAKE_shadow(splitstr_cmd,var_map):
    assert splitstr_cmd[1] == "shadow" 
    assert splitstr_cmd[2] == "with" 

    q = splitstr_cmd[3].split(",")
    assert len(q) == 3 

    assert q[0] in var_map 
    prg = var_map[q[0]] 
    prg = MAIN_method_for_object(prg) 
    file_path = q[2] 
    fitting_struct = q[1] 
    assert os.path.exists(q[2])

    sg = ShadowGen(prg,file_path,fitting_struct)
    return sg 

def MAKE_ssino(splitstr_cmd,var_map):
    assert splitstr_cmd[1] == "ssino" 
    assert splitstr_cmd[2] == "with" 

    q = splitstr_cmd[3].split(",")
    l = len(q)
    assert l in {3,4,5,6,7}
    
    r0 = int(q[0])
    
    assert q[1] in var_map 
    prg = var_map[q[1]]
    prg = MAIN_method_for_object(prg) 

    assert q[2] in var_map 
    prg2 = var_map[q[2]] 
    prg2 = MAIN_method_for_object(prg2) 

    lcg_input_only = 0 
    uniform_io_dist = 1 
    shuffle_dist = 0 
    prg_io_noise = 1 

    if l >= 4: 
        lcg_input_only = int(q[3])
        assert lcg_input_only in {0,1}

    if l >= 5: 
        uniform_io_dist = int(q[4])
        assert uniform_io_dist in {0,1} 

    if l >= 6: 
        shuffle_dist = int(q[5])
        assert shuffle_dist in {0,1} 

    if l >= 7: 
        prg_io_noise = int(q[6])
        assert prg_io_noise in {0,1}  
    
    sno = SSINetOp.one_instance(r0,prg,prg2,\
        lcg_input_only,uniform_io_dist,shuffle_dist,prg_io_noise)

    return sno  

# TODO: incomplete 
def MAKE_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make"

    if "lcg" in splitstr_cmd[1]: 
        return MAKE_lcgvx(splitstr_cmd,var_map)

    if "mdr" in splitstr_cmd[1]: 
        return MAKE_mdrvx(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "optri":
        assert splitstr_cmd[2] == "with" 
        
        parameters = splitstr_cmd[3].split(",") 
        assert len(parameters) == 5 

        int_seed = int(parameters[0])

        assert parameters[1] in var_map 
        prg = var_map[parameters[1]] 
        prg = MAIN_method_for_object(prg) 
        
        prg_ = prg__single_to_int(prg)
        gen_type = int(parameters[2]) 
        assert gen_type in {1,2} 

        add_noise = bool(int(parameters[3]))

        assert parameters[4] in var_map
        base_sequence = var_map[parameters[4]]
        assert type(base_sequence) == list and len(base_sequence) >= 2 
        base_sequence = IntSeq(base_sequence) 
        M = base_sequence.difftri(cast_type=np.int32)

        return OpTriGen(int_seed,M,prg_,gen_type,forward_func=add,\
            backward_func=sub,add_noise=add_noise) 

    if splitstr_cmd[1] == "multimetric": 
        assert splitstr_cmd[2] == "with" 
        assert splitstr_cmd[3] in var_map
        lx = var_map[splitstr_cmd[3]] 
        return MultiMetric(lx)

    if splitstr_cmd[1] == "qval":
        assert splitstr_cmd[2] == "with" 

        parameters = splitstr_cmd[3].split(",")

        assert parameters[0] in var_map 
        V = var_map[parameters[0]] 

        assert parameters[1] in var_map 
        index_selector = var_map[parameters[1]] 

        assert parameters[2] in var_map 

        def length_outputter(): 
            q = var_map[parameters[2]] 
            return int(round(q()))

        assert parameters[3] in var_map 
        range_outputter = var_map[parameters[3]] 

        try: 
            adj_type = int(parameters[4])
        except: 
            raise ValueError("invalid adjustment type {}".format(parameters[4])) 

        assert adj_type in {1,2} 

        return QValueOutputter(IntSeq(V),index_selector,length_outputter,\
            range_outputter,adj_type) 

    if splitstr_cmd[1] == "pid":
        assert splitstr_cmd[2] == "with" 

        parameters = splitstr_cmd[3].split(",")
        assert parameters[0] in var_map 
        prg = var_map[parameters[0]] 
        prg = MAIN_method_for_object(prg) 
    
        assert parameters[1] in var_map 
        f_out = var_map[parameters[1]] 
        f_out = MAIN_method_for_object(f_out) 
        def f_out_(): return int(round(f_out()))

        assert parameters[2] in var_map 
        l_out = var_map[parameters[2]] 
        l_out = MAIN_method_for_object(l_out) 
        def l_out_(): return int(round(l_out()))

        assert parameters[3] in var_map 
        r_out = var_map[parameters[3]] 
        r_out = MAIN_method_for_object(r_out)

        adjustment_type = int(parameters[4]) 
        assert adjustment_type in {1,2} 

        return PIDValueOutputter(prg,f_out_,\
            l_out_,r_out,adjustment_type)

    if splitstr_cmd[1] == "op2": 
        return MAKE_op2(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "rch":
        return MAKE_rch(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "iomaps":
        return MAKE_iomaps(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "echo":
        return MAKE_echo(splitstr_cmd) 

    if splitstr_cmd[1] == "shadow":
        return MAKE_shadow(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "idforest": 
        return MAKE_idforest(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "ssino": 
        return MAKE_ssino(splitstr_cmd,var_map) 

    raise ValueError("diffektor") 
