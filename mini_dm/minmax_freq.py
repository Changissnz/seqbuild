from collections import defaultdict 

"""
given a dictionary with string-like key and value as list, <MinMaxFreq> 
calculates the frequency of each element in the value. 
"""
class MinMaxFreq: 

    def __init__(self,d:dict,log_revd:bool): 
        self.d = d  
        self.c = defaultdict(int)
        self.revd = defaultdict(list)
        self.log_revd = log_revd
        self.sorted_counts = None
        self.current_key = None 
        self.fin = False 
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
        self.fin = False 
        self.add_pair_update()  

    def add_pair_update(self,k):
        v = self.d[k] 
        for v_ in v: 
            self.c[v_] += 1 
            if self.log_revd: 
                self.revd[v_].append(k) 
        v.clear() 
        return

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
        self.sorted_counts = sorted(self.c,key=lambda x:x[1])
        return