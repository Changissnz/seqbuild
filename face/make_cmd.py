from seqgen.gg_gen import * 
from intigers.lcg_v3 import * 

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","idforest","optri",\
             "rch","qval"] 

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
            exclude_zero__auto_td = bool(parameters[8])
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

# TODO: incomplete 
def MAKE_proc(splitstr_cmd,var_map): 
    assert splitstr_cmd[0] == "make"

    if "lcg" in splitstr_cmd[1]: 
        return MAKE_lcgvx(splitstr_cmd,var_map)

    if splitstr_cmd[1] == "mdr": 
        assert splitstr_cmd[2] == "with" 
        assert splitstr_cmd[3] in var_map 
        lx = var_map[splitstr_cmd[3]] 
        mdx = ModuloDecomp(IntSeq(lx)) 
        mdx.merge(False)
        return ModuloDecompRepr(mdx,reconstruct_type=1)

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
        length_outputter = var_map[parameters[2]] 

        assert parameters[3] in var_map 
        range_outputter = var_map[parameters[3]] 

        try: 
            adj_type = int(parameters[4])
        except: 
            raise ValueError("invalid adjustment type {}".format(parameters[4])) 

        assert adj_type in {1,2} 

        return QValueOutputter(V,index_selector,length_outputter,\
            range_outputter,adj_type) 

    return None
