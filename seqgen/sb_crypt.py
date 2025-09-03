from .rch_gen import * 
from types import MethodType,FunctionType
import os 

DEFAULT_CHARACTER_MAP = ord 
DEFAULT_CHARACTER_IMAP = chr 

"""
?one-way encryption? 
"""
class SBCrypt:

    def __init__(self,filepath,out_filepath,prg,uncertainty_upper_threshold,\
        char_mapper=DEFAULT_CHARACTER_MAP,ichar_mapper=DEFAULT_CHARACTER_IMAP):  
        assert os.path.exists(filepath)
        assert type(prg) in {MethodType,FunctionType}

        self.file_path = filepath
        self.out_filepath = out_filepath
        self.fi_obj,self.fi_end = None,None 
        self.set_fileobj()
        self.prg = prg 
        self.uut = uncertainty_upper_threshold

        self.cmap = char_mapper 
        self.icmap = ichar_mapper 

        self.invalid_char_stat = False 
        self.fin_stat = False 

        self.charseq = [] 
        self.c_index = 0 
        return

    def encrypt_one_chunk(self): 
        m = int(round(modulo_in_range(self.prg(),[3,50])))

        self.read_n_characters(m) 
        if len(self.charseq) == 0: 
            self.fin_stat = True 
            return 

        cseq = self.map_charseq()
        for c in cseq: 
            if self.c_index >= 10: 
                self.c_index = 0 
                self.fi_obj2.write("\n") 
            self.fi_obj2.write(str(c) + ",")
            self.c_index += 1 
        return

    def read_n_characters(self,n): 
        self.charseq.clear()
        for _ in range(n):
            if self.fi_obj.tell() == self.fi_end: 
                break 
            char = self.fi_obj.read(1)
            self.charseq.append(char)
        return

    def map_charseq(self):

        prg_ = prg__single_to_int(self.prg)
        num_nodes = 4
        dim_range = [3,5] 
        ufreq_rate = [3,9]
        rchv = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg_,\
                ufreq_rate,mutrate=0.5,queue_capacity=1000)
        rch = rchv.rch

        cseq = [] 
        for c in self.charseq: 
            q,qx = None,None 
            try: 
                q = self.map_char(c) 
                qx = rch.apply(q) 
            except: 
                self.invalid_char_stat = True 
            if type(qx) != type(None): 
                qx = modulo_in_range(qx,(-5000,5000.1))
            cseq.append(qx)  
        return cseq 

    def map_char(self,c):
        q = None 
        try:
            q = self.cmap(c) 
        except:
            print("no can do.") 
            self.invalid_char_stat = True 
            return 
        return q  

    def set_fileobj(self):
        self.fi_obj = open(self.file_path,'r')
        self.fi_obj.seek(0,os.SEEK_END) 
        self.fi_end = self.fi_obj.tell()
        self.fi_obj.seek(0)
        self.fi_obj2 = open(self.out_filepath,'w')
        return 

    def close(self): 
        self.fi_obj.close() 
        self.fi_obj2.flush() 
        self.fi_obj2.close() 
        
        del self 

