from desi.multi_metric import * 

class AGV2GuidedGen: 

    def __init__(self,base_prg,aux_prg):
        self.base_prg = base_prg 
        self.aux_prg = aux_prg 

        self.spannistia = None 
        self.output_mode = "base" 
        return 

    def setspan(self,span):
        self.spannistia = span 