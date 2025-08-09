"""
base file for interface language design and corresponding 
interpretive functions to process commands. 
"""
from .cmd_proc import * 

class CommLangParser:

    def __init__(self,filepath:str):
        self.fp = filepath
        self.cmdlines = []
        self.vartable = dict() 

    def load_next_lines(self): 
        return -1 

    def load_cmd_lines(self,cmdlines):
        return -1 

    def SET_proc(self,splitstr_cmd):
        assert splitstr_cmd[0] == "set" 
        n = splitstr_cmd[1] 
        assert splitstr_cmd[2] == "=" 

        # make, run , open

        # case 
        if splitstr_cmd[3] == "make":
            self.vartable[n] = MAKE_proc(splitstr_cmd[3:])
        elif splitstr_cmd[3] == "run": 
            self.vartable[n] = RUN_proc(splitstr_cmd[3:],self.vartable) 
        elif splitstr_cmd[3] == "open": 
            self.vartable[n] = OPEN_proc(splitstr_cmd[3:])
            return
        else: 
            assert False 

        return n, self.vartable[n]