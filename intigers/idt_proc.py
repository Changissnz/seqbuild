from .idectree import * 
from types import MethodType 

class IDTProc:

    def __init__(self,tn,proc_mode="T"): 
        assert type(tn) == IDecNode
        assert proc_mode in {"T","O+T"}
        self.tn = tn
        self.proc_mode = proc_mode
        self.flow_queue = set() 

    #----------------------- input-flow processing functions 

    def __next__(self):
        nfq = set() 
        vx = [] 
        dxo = defaultdict(set) 
        opu = lambda x,x2: x | x2 

        while len(self.flow_queue) > 0:
            tn = self.flow_queue.pop()
            vs,dx = self.process_node_queue(tn)

            vx.extend(vs) 
            dxo = numberdict_op(dxo,dx,opu)

            for k in dx.keys():
                tn2 = tn.fetch_conn(k)
                nfq = nfq | {tn2}

        self.flow_queue = nfq
        self.update_node_queue(dxo)
        return

    def process_node_queue(self,tn):
        assert type(tn) == IDecNode

        vs = []
        dx = defaultdict(set) 
        for x in tn.acc_queue:
            q = tn.entryf(x) 
            vs.append(q)

            q2 = tn.travf.apply(x) 
            if type(q2) != type(None):
                dx[q2] |= {x}

        tn.clear() 
        return vs,dx

    def update_node_queue(self,dx):

        for k,v in dx.items():
            tn = None
            for x in self.flow_queue:
                if x.idn == k:
                    tn = x
                    break
            assert type(tn) != type(None) 
            tn.add_to_acc_queue(list(v))
        return

    def inflow_set(self,iset):
        assert type(iset) == iset
        self.tn.add_to_acc_queue(list(iset))

        if self.tn not in self.flow_queue:
            self.flow_queue |= {self.tn}
        return

    #----------------------- standalone value-processing functions

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