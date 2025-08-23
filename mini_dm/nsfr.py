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

    def __init__(self,file_obj,cast_func,is_periodic:bool=True): 
        assert type(file_obj) == io.TextIOWrapper

        self.fi_obj = file_obj 
        self.cast_func = cast_func 
        self.is_periodic = is_periodic
        self.is_finished = False 
        self.set_fileobj()

        self.cast_func = cast_func 
        self.cache = []

    def __next__(self):
        if len(self.cache) > 0: 
            return self.cache.pop(0)

        stat = True 
        while stat: 
            stat = self.read_one_line()
            if len(self.cache) > 0: break 
        
        if len(self.cache) == 0: 
            return None 
        return self.cache.pop(0)

    def read_one_line(self): 
        if self.is_finished: 
            return False 

        if self.fi_obj.tell() == self.fi_end: 
            self.is_finished = True 
            self.reload() 

        l = self.fi_obj.readline().strip()
        if len(l) == 0: 
            return True 

        V = string_to_vector(l,self.cast_func)
        self.cache.extend(V) 
        return True 

    def reload(self):
        if not self.is_periodic:
            return 
        self.set_fileobj()
        self.is_finished = False 
        return

    def set_fileobj(self):
        self.fi_obj.seek(0,os.SEEK_END) 
        self.fi_end = self.fi_obj.tell()
        self.fi_obj.seek(0)
        return 

    def close(self): 
        self.fi_obj.close() 
        del self 