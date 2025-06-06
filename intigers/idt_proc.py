from .idectree import * 
from types import MethodType 

class IDTProc:

    def __init__(self,tn,proc_mode="T"): 
        assert type(tn) == IDecNode
        assert proc_mode in {"T","O+T"}
        self.tn = tn
        self.proc_mode = proc_mode

    """
    produce a path that v undergoes starting at `tn`. 
    """
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

    """
    sequence of output values from the `entryf` functions applied 
    on `v` if `v` were to travel the sequence of nodes `vpath`. 
    """
    def travel_path(self,v,vpath):

        xr = self.tn 
        def add_entry(v_): 
            # print out first value 
            vx = xr.entryf(v_)
            if type(vx) == np.ndarray:
                vx = vx.flatten()
            else: 
                vx = [vx]
            V.extend(vx)

        V = []   
        add_entry(v) 

        for v2 in vpath: 
            xr = xr.fetch_conn(v2)
            add_entry(v2)
        return V

    """
    travels `v2` on the path `p` that `v1` 
    takes from `process_value`. Outputs the
    corresponding `entryf` values, path-isomorphic 
    to the output of `v1`. 
    """
    def iso_output(self,v1,v2):
        p = self.process_value(v1)
        return self.travel_path(v2,p[1:])

    """
    outputs the map `node idn` -> `number of input elements that traveled the node`
    """
    def traffic_map(self,S):
        D = defaultdict(int)
        for s in S:
            pth = self.process_value(v)
            pth = dict({(p,1) for p in pth})
            D = numberdict_op(D,pth,add)
        return D  