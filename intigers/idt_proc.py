from .idectree import * 

class IDTProc:

    def __init__(self,tn,proc_mode="T"): 
        assert type(tn) == IDecNode
        assert proc_mode in {"T","O+T"}
        self.tn = tn
        self.proc_mode = proc_mode

    def process_value(self,v):
        tn = self.tn 
        path = [tn.idn] 
        while True: 
            if type(tn.travf) == type(None): 
                break 
            q = tn.travf.apply(v)
            if type(q) == type(None): break
            tn2 = tn.fetch_conn(q)
            path.append(tn2.idn)
            tn = tn2 
        return path 
