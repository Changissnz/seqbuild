from seqgen.gg_gen import * 
from seqgen.mdr_gen import * 
from seqgen.optri_gen import * 
from seqgen.rch_gen import * 
from seqgen.idt_gen import * 
from seqgen.shadow_gen import * 
from seqgen.n2mv_gen import * 
from seqgen.gg_gen import * 
from seqgen.afs_gen import * 
from seqgen.fit22_gen import * 
from seqgen.lps_gen import * 
from desi.fraction import * 
from desi.differentials import * 
from desi.prng_dec_delta import * 
from intigers.lcg_v3 import * 
from intigers.prng_pw_op import *
from seqgen.sb_crypt import * 
from seqgen.ssi_netop import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","mdrv2",\
        "mdrgen","iomaps","idforest","optri","rch","qval",\
        "pid","echo","shadow","ssino","n2m","gg","afs","fit22",\
        "pddelta"]   

ARITHMETIC_OP_STR_MAP = {"+":add,"-":sub,"/":zero_div0,"*":mul,"%":safe_mod}  

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

    if type(q) == PRNGDecimalDelta: 
        return q.__next__ 

    if type(q) == N2MVGen: 
        return q.__next__

    if type(q) == GaugeGuidedGen: 
        return q.__next__ 

    if type(q) == AFSGen: 
        return q.__next__

    if type(q) == Fit22Gen: 
        return q.__next__ 

    if type(q) == LPSGen:
        return q.__next__ 

    return -1 
    assert False, "{} is not valid".format(q) 


#------------------------------------------- class 3 PRNGs: stacked vector-based form-fitters + noise adders 

"""
prg
primary generator 

prg,prg
primary generator, secondary generator 

prg,prg,range
primary generator, secondary generator, NM-dimension range

prg,prg,range,range
primary generator, secondary generator, NM-dimension range, update frequency range 
"""
def MAKE_n2m(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make" 
    assert "n2m" == splitstr_cmd[1] 
    assert "with" == splitstr_cmd[2] 

    q = splitstr_cmd[3].split(",")
    l = len(q)
    assert 1 <= l <= 6, "invalid parameters"

    assert q[0] in var_map 
    prg = var_map[q[0]] 
    prg = MAIN_method_for_object(prg) 

    if l == 1: 
        return N2MVGen(DEFAULT_N2MV_NM_DIMENSION_RANGE,prg,None,DEFAULT_N2MV_NM_UPDATE_FREQUENCY_RANGE)

    assert q[1] in var_map 
    prg2 = MAIN_method_for_object(var_map[q[1]]) 
    if l == 2: 
        return N2MVGen(DEFAULT_N2MV_NM_DIMENSION_RANGE,prg,prg2,DEFAULT_N2MV_NM_UPDATE_FREQUENCY_RANGE)

    assert l in {4,6}, "invalid parameters" 

    r0,r1 = int(q[2]),int(q[3])
    assert is_valid_range((r0,r1),True,False) 

    if l == 4: 
        return N2MVGen((r0,r1),prg,prg2,DEFAULT_N2MV_NM_UPDATE_FREQUENCY_RANGE)

    s0,s1 = int(q[4]),int(q[5])
    return N2MVGen((r0,r1),prg,prg2,(s0,s1))

"""
prg, int
primary generator, vector size 

prg, prg2, int
primary generator, secondary generator, vector size 

prg, prg2, prg3, int
primary generator, secondary generator, third generator, vector size 

prg, prg2, prg3, prg4, int
primary generator, secondary generator, third generator, fourth generator, vector size 

"""
def MAKE_afs(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make" 
    assert "afs" == splitstr_cmd[1] 
    assert "with" == splitstr_cmd[2] 

    q = splitstr_cmd[3].split(",")
    l = len(q)
    assert 2 <= l <= 5, "invalid parameters"

    prg = MAIN_method_for_object(var_map[q[0]]) 

    if l == 2: 
        x = int(q[1]) 
        return AFSGen(prg,None,None,None,x,True) 

    prg2 = MAIN_method_for_object(var_map[q[1]]) 

    if l == 3: 
        x = int(q[2]) 
        return AFSGen(prg,prg2,None,None,x,True) 

    prg3 = MAIN_method_for_object(var_map[q[2]]) 

    if l == 4: 
        x = int(q[3]) 
        return AFSGen(prg,prg2,prg3,None,x,True) 

    prg4 = MAIN_method_for_object(var_map[q[3]]) 
    x = int(q[4]) 
    return AFSGen(prg,prg2,prg3,prg4,x,True) 

#------------------------------------------- class D PRNGs: differentials of fitters between pairs of two-dimensional points 


def MAKE_fit22(splitstr_cmd,var_map): 

    assert splitstr_cmd[0] == "make" 
    assert "fit22" == splitstr_cmd[1] 
    assert "with" == splitstr_cmd[2] 

    q = splitstr_cmd[3].split(",")
    l = len(q)
    assert 1 <= l <= 4, "invalid parameters"

    prg = MAIN_method_for_object(var_map[q[0]]) 

    if l == 1: 
        return Fit22Gen(prg,None,None,None) 

    prg2 = MAIN_method_for_object(var_map[q[1]]) 

    if l == 2: 
        return Fit22Gen(prg,prg2,None,None) 

    prg3 = MAIN_method_for_object(var_map[q[2]]) 

    if l == 3: 
        return Fit22Gen(prg,prg2,prg3,None) 

    prg4 = MAIN_method_for_object(var_map[q[3]]) 
    return Fit22Gen(prg,prg2,prg3,prg4)

def MAKE_lps(splitstr_cmd,var_map): 

    assert splitstr_cmd[0] == "make" 
    assert "lps" == splitstr_cmd[1] 
    assert "with" == splitstr_cmd[2]

    q = splitstr_cmd[3].split(",")
    l = len(q)
    assert 1 <= l <= 3, "invalid parameters"

    prg = MAIN_method_for_object(var_map[q[0]]) 
    if l == 1: 
        return LPSGen(prg,None,None) 

    prg2 = MAIN_method_for_object(var_map[q[1]]) 
    if l == 2: 
        return LPSGen(prg,prg2,None)

    prg3 = MAIN_method_for_object(var_map[q[2]]) 
    return LPSGen(prg,prg2,prg3) 

#-------------------------------------------- class 0 PRNGs: LCGs and LCG form-fitters 

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

        assert len(parameters) in {6,7,9,11,13},"parameters len {}: {}".format(len(parameters),parameters)

        px = [float(p) for p in parameters[:5]] 
        prg = MAIN_method_for_object(var_map[parameters[5]]) 

        if len(parameters) == 6:
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],prg,super_range=None)
        elif len(parameters) == 7: 
            add_noise = bool(int(parameters[6])) 
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],prg,super_range=None,add_noise=add_noise)

        elif len(parameters) == 9: 
            r0,r1 = float(parameters[6]),float(parameters[7])
            add_noise = bool(int(parameters[8])) 
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],prg,super_range=(r0,r1),add_noise=add_noise)

        elif len(parameters) == 11: 
            r0,r1 = float(parameters[6]),float(parameters[7])
            ternary_size_range = (int(parameters[8]),int(parameters[9]))
            add_noise = bool(int(parameters[10]))
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],prg,super_range=(r0,r1),ternary_size_range=ternary_size_range,\
                add_noise=add_noise)
        else: 
            r0,r1 = float(parameters[6]),float(parameters[7])
            ternary_size_range = (int(parameters[8]),int(parameters[9]))
            ternary_delta_timestamp_range = (int(parameters[8]),int(parameters[9]))

            add_noise = bool(int(parameters[12]))
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],prg,super_range=(r0,r1),ternary_size_range=ternary_size_range,\
                ternary_delta_timestamp_range = ternary_delta_timestamp_range,add_noise=add_noise)
        return g3
        
    assert False 

def MAKE_mdrvx(splitstr_cmd,var_map): 
    assert splitstr_cmd[2] == "with" 
    parameters = splitstr_cmd[3] 
    parameters = parameters.split(",") 
    assert parameters[0] in var_map 

    if splitstr_cmd[1] == "mdr": 
        assert 1 <= len(parameters) <= 2
        lx = var_map[parameters[0]] 
        max_absmult = None 

        if len(parameters) == 2: 
            max_absmult = int(parameters[1]) 
        mdx = ModuloDecomp(IntSeq(lx),max_absmult)  
        mdx.merge(False)
        return ModuloDecompRepr(mdx,reconstruct_type=1)

    if splitstr_cmd[1] == "mdrv2":
        assert len(parameters) < 4

        exclude_neg = False 
        max_absmult = None 

        if len(parameters) >= 2:
            i = bool(int(parameters[1]))
            exclude_neg = i 
        
        if len(parameters) == 3: 
            max_absmult = int(parameters[2])     
        
        lx = var_map[parameters[0]]
        mdx = ModuloDecompV2(IntSeq(lx),exclude_neg,max_absmult)
        return ModuloDecompRepr(mdx,reconstruct_type=2) 

    if splitstr_cmd[1] == "mdrgen": 
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

#-------------------------------------------- class 1 PRNGs: linear combination/polynomial based 
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
    output_range = [-10000,10000]
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

    idf = IDecForest(IntSeq(V),MO,cache_size,reprod_rate_range,max_trees,G,prg2,True,False)   
    return idf

#-------------------------------------------------------- class I: identity PRNGs (cycling of file)

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

#---------------------------------------------------------- class 2: multi-factored form-fitters + noise-adders 

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

"""
prg,range 
primary generator, super-range 

prg,range, 0|1
primary generator, super-range, mode::(cov|uwpd) 

prg,range, 0|1,0|1
primary generator, super-range, mode::(cov|uwpd),allow subrange drift

prg,range, 0|1,0|1,range
primary generator, super-range, mode::(cov|uwpd),allow subrange drift, base output span 

prg,range, 0|1,0|1,range,density measure 
primary generator, super-range, mode::(cov|uwpd),allow subrange drift, base output span 
"""
def MAKE_gg(splitstr_cmd,var_map): 
    assert splitstr_cmd[1] == "gg" 
    assert splitstr_cmd[2] == "with" 

    q = splitstr_cmd[3].split(",")
    assert len(q) >= 3 

    assert q[0] in var_map 
    prg = MAIN_method_for_object(var_map[q[0]]) 
    r0,r1 = int(float(q[1])),int(float(q[2])) 

    r0 = modulo_in_range(r0,DEFAULT_AGV2GG_BASE_OUTPUT_SPAN)
    r1 = modulo_in_range(r1,DEFAULT_AGV2GG_BASE_OUTPUT_SPAN)
    if r0 == r1: r1 += 1 
    r0,r1 = sorted([r0,r1])
    if len(q) == 3: 
        return GaugeGuidedGen(prg,DEFAULT_AGV2GG_BASE_OUTPUT_SPAN,DEFAULT_AGV2GG_DENSITY_MEASURE,\
            (r0,r1),target_measure="cov",allow_subrange_drift=True) 

    x =  int(q[3]) 
    assert x in {0,1}
    target_measure = "cov" if x == 0 else "uwpd"

    if len(q) == 4: 
        return GaugeGuidedGen(prg,DEFAULT_AGV2GG_BASE_OUTPUT_SPAN,DEFAULT_AGV2GG_DENSITY_MEASURE,\
            (r0,r1),target_measure=target_measure,allow_subrange_drift=True) 

    x = int(q[4]) 
    assert x in {0,1}
    allow_subrange_drift = bool(int(x)) 

    if len(q) == 5: 
        return GaugeGuidedGen(prg,DEFAULT_AGV2GG_BASE_OUTPUT_SPAN,DEFAULT_AGV2GG_DENSITY_MEASURE,\
            (r0,r1),target_measure=target_measure,allow_subrange_drift=allow_subrange_drift) 

    s0,s1 = int(q[5]),int(q[6]) 
    if len(q) == 7: 
        return GaugeGuidedGen(prg,(s0,s1),DEFAULT_AGV2GG_DENSITY_MEASURE,\
            (r0,r1),target_measure=target_measure,allow_subrange_drift=allow_subrange_drift) 

    density_measure = int(q[7]) 

    if len(q) == 8: 
        return GaugeGuidedGen(prg,(s0,s1),density_measure,\
            (r0,r1),target_measure=target_measure,allow_subrange_drift=allow_subrange_drift) 
    assert False, "invalid parameters" 

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

#---------------------------------------------- class T PRNGs: time-based

def MAKE_pddelta(splitstr_cmd,var_map): 
    assert splitstr_cmd[1] == "pddelta" 
    assert splitstr_cmd[2] == "with" 

    r1 = float(splitstr_cmd[3]) 
    return PRNGDecimalDelta(r1)  


#--------------------------------------------- extra structures for use as PRNG parameters 


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

#---------------------------------------------- class V PRNGs: quotient and polynomial form-fitters of vector inputs. 

def MAKE_qval(splitstr_cmd,var_map): 
    assert splitstr_cmd[1] == "qval" 
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

def MAKE_pid(splitstr_cmd,var_map): 

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

#--------------------------------------------- class O PRNGs: pairwise arithmetic operations on vector inputs. 

def MAKE_optri(splitstr_cmd,var_map): 
    
    assert splitstr_cmd[1] == "optri" 
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

#------------------------------------------------------------------------------------- 

# TODO: incomplete 
def MAKE_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make"

    if "lcg" in splitstr_cmd[1]: 
        return MAKE_lcgvx(splitstr_cmd,var_map)

    if "mdr" in splitstr_cmd[1]: 
        return MAKE_mdrvx(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "optri":
        return MAKE_optri(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "multimetric": 
        assert splitstr_cmd[2] == "with" 
        assert splitstr_cmd[3] in var_map
        lx = var_map[splitstr_cmd[3]] 
        return MultiMetric(lx)

    if splitstr_cmd[1] == "qval":
        return MAKE_qval(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "pid":
        return MAKE_pid(splitstr_cmd,var_map)

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

    if splitstr_cmd[1] == "n2m": 
        return MAKE_n2m(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "gg": 
        return MAKE_gg(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "afs": 
        return MAKE_afs(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "fit22":
        return MAKE_fit22(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "lps": 
        return MAKE_lps(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "pddelta": 
        return MAKE_pddelta(splitstr_cmd,var_map) 

    raise ValueError("diffektor") 
