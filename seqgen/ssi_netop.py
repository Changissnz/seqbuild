from .ssi_load import * 
from morebs2.numerical_generator import prg_partition_for_sz__n_rounds

class SSINetOp:

    def __init__(self,struct_list,h2tree_map,prg):
        assert len(struct_list) > 3 
        assert len(h2tree_map) > 0
        assert type(prg) in {MethodType,FunctionType}

        self.struct_list = struct_list
        self.h2tree_map = h2tree_map
        self.prg = prg 

        # storage of values from some source; values can be used for 
        # the structures `mdr`,`mdrv2`,`optri`. 
        self.mainstream_queue = []
        self.tree_idns = sorted(self.h2tree_map.keys())
        self.t_index = 0  
        return 

    @staticmethod 
    def one_instance(num_nodes,prg,prg2): 
        assert num_nodes >= 5 

        num_sets = 3 
        var = 0.25 
        num_rounds = 10 
        prt = prg_partition_for_sz__n_rounds(num_nodes,\
            num_sets,prg,var,num_rounds)
        slist = ["lcg","mdr","optri"] 
        
        sls = []

        # make each of the types 
        for (s,p) in zip(slist,prt): 
            sl = SSINetOp.one_instance__slist(s,p,prg)
            sls.extend(sl)

        lx = len(sls)
        q = max([round(lx/3),5])
        q = int(q)
        ssb2n = SSIBatch2Net(sls,q,prg) 
        ssb2n.make_net()

        q1 = sls 
        q2 = ssb2n.h2tree_map
        return SSINetOp(q1,q2,prg2)

    @staticmethod
    def one_instance__slist(sidn,batch_size,prg):
        # make an LCG from `prg` output 
        ix = [prg() for _ in range(4)] 
        ix = [modulo_in_range(ix_,[1.,101.]) for ix_ in ix]
        prg2 = prg__LCG(ix[0],ix[1],ix[2],ix[3])

        if sidn == "lcg":
            prg3 = prg__single_to_bounds_outputter(prg,4) 
            param_bounds = prg3()
        elif sidn == "mdr": 
            param_bounds = [0,1]
        else:
            param_bounds = None 

        ssibl = SSIBatchLoader__TypeLCGNet(sidn,param_bounds,batch_size,prg)
        ssibl.instantiate_slist()
        return ssibl.slist 
