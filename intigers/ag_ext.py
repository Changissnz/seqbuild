from morebs2.aprng_gauge import *

class APRNGGaugeV2(APRNGGauge):

    def __init__(self,aprng,frange,pradius:float):
        super().__init__(aprng,frang,pradius)

    def measure_cycle(self,max_size,\
        term_func=lambda l,l2: type(l) == type(None)): 
        
        super().measure_cycle(max_size,term_func)