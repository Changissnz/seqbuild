from .ssi_load import * 


class SSINetOp:

    def __init__(self,struct_list,h2tree_map):
        self.struct_list = struct_list
        self.h2tree_map = h2tree_map

        # storage of values from some source; values can be used for 
        # the structures `mdr`,`mdrv2`,`optri`. 
        self.mainstream_queue = []
        return 