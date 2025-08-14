from seqgen.gg_gen import * 
from seqgen.mdr_gen import * 
from seqgen.optri_gen import * 
from desi.fraction import * 
from desi.differentials import * 
from intigers.lcg_v3 import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","mdrv2",\
        "mdrgen","idforest","optri","rch","qval",\
        "pid"]  


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
            q.reset_first(int(first))
            return q.reconstruct() 
        return f 

    if type(q) == MDRGen: 
        return q.__next__ 

    if type(q) == OpTriGen: 
        return q.__next__ 

    if type(q) == PIDValueOutputter: 
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

            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],0,False,prg,super_range,\
                exclude_zero__auto_td,is_rmod)

        # case: trinary guided w/ option to exclude 0 from trinary vector guides. 
        elif len(parameters) == 9: 
            exclude_zero__auto_td = bool(int(parameters[8]))
            is_rmod = False 

            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],0,False,prg,super_range,\
                exclude_zero__auto_td,is_rmod)

        # case: trinary guided w/ option to exclude 0 from trinary vector guides, 
        # option to use `reflective modification` output mode 
        else: 
            exclude_zero__auto_td = bool(int(parameters[8]))
            is_rmod = bool(int(parameters[9]))
            g3 = LCGV3(px[0],px[1],px[2],px[3],px[4],0,False,prg,super_range,\
                exclude_zero__auto_td,is_rmod)
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

# TODO: incomplete 
def MAKE_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make"

    if "lcg" in splitstr_cmd[1]: 
        return MAKE_lcgvx(splitstr_cmd,var_map)

    if "mdr" in splitstr_cmd[1]: 
        return MAKE_mdrvx(splitstr_cmd,var_map) 

    if splitstr_cmd[1] == "optri":
        '''
        int seed,prg,gen type, add noise,base sequence (for m)
        '''
        assert splitstr_cmd[2] == "with" 
        
        parameters = splitstr_cmd[3].split(",") 
        assert len(parameters) == 5 

        int_seed = int(parameters[0])

        assert parameters[1] in var_map 
        prg = var_map[parameters[1]] 
        prg = MAIN_method_for_object(prg) 
        
        def prg_(): return int(round(prg())) 
        
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

    return None
