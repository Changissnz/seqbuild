import io,os 
from morebs2.matrix_methods import string_to_vector

"""
(N)umerical (S)equence (F)ile (R)eader. 

Reads a comma-separated value (csv) file object `file_obj` to 
load values, usually integers or floats (as determined by the 
`castFunc` argument), line by line into the output `cache`. 

File reader effectively serves as a numerical generator by reading 
from a file that contains a numerical sequence, separated by commas 
and newlines. 
"""
class NSFileReader: 

    def __init__(self,file_obj,castFunc): 
        assert type(file_obj) == io.TextIOWrapper
        
        self.f_obj = file_obj 
        self.f_obj.seek(0,os.SEEK_END) 
        self.f_end = self.f_obj.tell()
        self.f_obj.seek(0)  

        self.castFunc = castFunc 
        self.cache = []

        self.fin_stat = False 

    def __next__(self):
        if self.fin_stat: return None 

        if len(self.cache) == 0:
            self.load_next_sequence() 
            return self.__next__() 
        q = self.cache.pop(0)
        return q 

    def load_next_sequence(self):
        l,stat = self.load_next_line() 

        if not stat: 
            self.fin_stat = True
            return 

        svec = string_to_vector(l,self.castFunc) 
        self.cache.extend(svec) 

    def load_next_line(self):
        if self.f_obj.tell() != self.f_end:
            line = self.f_obj.readline()
            return line,True 
        return None,False 

    def close(self): 
        self.f_obj.close() 
        del self 