from .idectree import * 
from types import MethodType 

class IDTProc:

    def __init__(self,tn,proc_mode="T"): 
        assert type(tn) == IDecNode
        assert proc_mode in {"T","O+T"}
        tn.full_clear() 
        self.tn = tn
        self.proc_mode = proc_mode
        self.flow_queue = set() 

    #----------------------- input-flow processing functions 

    """
    Outputs a dictionary of 
        node idn. -> list of values from `entryf` function 
    """
    def __next__(self):
        nfq = set() 
        vx = defaultdict(list)
        dxo = defaultdict(set) 

        def opu(x,x2):
            if x == 0:
                return x2
            elif x2 == 0: 
                return x
            return x | x2

        while len(self.flow_queue) > 0:
            tn = self.flow_queue.pop()
            vs,dx = self.process_node_queue(tn)

            vx[tn.idn] = vs
            dxo = numberdict_op(dxo,dx,opu)

            for k in dx.keys():
                tn2 = tn.fetch_conn(k)
                nfq = nfq | {tn2}

        self.flow_queue = nfq
        self.update_node_queue(dxo)
        return vx 

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
        assert type(iset) == set
        self.tn.add_to_acc_queue(list(iset))
        self.flow_queue |= {self.tn}
        return

    #----------------------- standalone value-processing functions

    # TODO: test 
    """
    D is a map of 
        node idn -> list of input values. 

    return: 
    - list, element is (input,output)
    """
    def splat_process(self,D):
        ks = sorted(D.keys())
        _,_,nodes = TNode.dfs(self.tn,False,False,True,None,set(ks))
        dnodes = {}
        for n in nodes: dnodes[n.idn] = n 

        L = []
        for k in ks:
            v = D[k]
            assert type(v) == list 

            n = dnodes[k]
            for v_ in v:
                q = n.travf.apply(v_)
                L.append((k,q))
        return L  

    """
    produce a path that v undergoes starting at `tn`. 
    """
    def process_value(self,v,calculate_vpath:bool=False):
        tn = self.tn 
        path = [tn.idn] 
        vpath = []
        while True: 
            if type(tn.travf) == type(None): 
                break 
            
            if calculate_vpath:
                q_ = tn.entryf(v)
                vpath.append(q_)

            q = tn.travf.apply(v)
            if type(q) == type(None): break
            tn2 = tn.fetch_conn(q)
            path.append(tn2.idn)
            tn = tn2 

        return path,vpath 

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
        p,_ = self.process_value(v1)
        return self.travel_path(v2,p[1:])

    """
    outputs the map `node idn` -> `number of input elements that traveled the node`
    """
    def traffic_map(self,S):
        D = defaultdict(int)
        for s in S:
            pth,_ = self.process_value(v)
            pth = dict({(p,1) for p in pth})
            D = numberdict_op(D,pth,add)
        return D  