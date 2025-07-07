from collections import defaultdict 
from morebs2.matrix_methods import is_vector

def vec_to_frequency_map(V):
    assert is_vector(V)
    dx = defaultdict(int)
    for v_ in V: dx[v_] += 1
    return dx 

def setseq_to_frequency_map(S): 
    dx = defaultdict(int)
    for s in S: 
        assert type(s) == set 
        for s_ in s: dx[s_] += 1 
    return dx 
"""
given a dictionary with string-like key and value as list, <MinMaxFreq> 
calculates the frequency of each element in the value. 
"""
class MinMaxFreq: 

    def __init__(self,d,log_revd:bool): 
        assert type(d) in {dict,defaultdict}
        if type(d) == dict: 
            d = defaultdict(list,d)
        self.d = d  
        self.c = defaultdict(int)
        self.revd = defaultdict(list)
        self.log_revd = log_revd
        self.sorted_counts = None
        self.current_key = None 
        self.fin = False 
        self.prev_fin = False 
        return

    def next_key(self):
        keyiter = iter(self.d.keys())
        k = None 

        while True: 
            try: 
                k = next(keyiter)
            except: 
                break

            if len(self.d[k]) > 0:
                return k 
            k = None
        self.current_key = k  
        return k  

    def add_pair(self,p,replace:bool=False): 
        assert len(p) == 2 and type(p[1]) == list 
        if not replace and p[0] in self.d: 
            print("[!] already present")
            return 
        self.d[p[0]] = p[1]
        self.prev_fin = self.fin 
        self.fin = False 
        self.add_pair_update()  

    def add_pair_update(self,k):
        v = self.d[k] 
        for v_ in v: 
            self.c[v_] += 1 
            if self.log_revd: 
                self.revd[v_].append(k) 
        v.clear() 
        self.fin = self.prev_fin 
        return

    """
    main method 
    """
    def count_one(self):
        if len(self.d[self.current_key]) == 0: 
            self.current_key = self.next_key() 
            if type(self.current_key) == type(None): 
                self.fin = True 
                return 
        element = self.d[self.current_key].pop(0) 
        self.c[element] += 1
        if self.log_revd: 
            self.revd[element].append(self.current_key) 

        for k,v in self.d.items(): 
            i = 0 
            while i < len(v): 
                if v[i] == element: 
                    self.c[element] += 1 
                    if self.log_revd: 
                        self.revd[element].append(k)
                    v.pop(i) 
                    continue 
                i += 1 

    def finalize_count(self):
        self.sorted_counts = sorted(self.c.items(),key=lambda x:x[1])
        self.fin = True 
        return

    # TODO: test 
    def nth_most_frequent(self,n):
        if type(self.sorted_counts) == type(None):
            return []
        if len(self.sorted_counts) == 0: 
            return [] 

        r = self.sorted_counts[-1][1]
        c = 0

        fx = []

        for x in self.sorted_counts[::-1]:
            if x[1] != r:
                c += 1
                r = x[1]
            if c > n: break

            if c == n:
                fx.append(x[0])  

        return fx 