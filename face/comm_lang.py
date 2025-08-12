"""
base file for interface language design and corresponding 
interpretive functions to process commands. 
"""
from .cmd_proc import * 

class CommLangParser:

    def __init__(self,filepath:str):
        self.fp = filepath
        self.file_obj = open(self.fp,"r") 

        self.file_obj.seek(0,os.SEEK_END) 
        self.file_end = self.file_obj.tell()
        self.file_obj.seek(0)  

        self.cmdlines = []
        self.commond = None 
        self.vartable = dict() 

    def close(self): 
        self.file_obj.close() 
        del self 

    def process_file(self): 
        stat = True

        while stat: 
            self.load_next_command()
            self.process_command() 
            stat = self.file_obj.tell() != self.file_end 

        return

    def load_next_command(self):
        comm = "" 
        stat = True 

        while stat: 
            if self.file_obj.tell() != self.file_end:
                line = self.file_obj.readline()
            else: 
                stat = False 
                continue 

            line = line.strip() 
            if len(line) == 0: 
                continue 

            # case: line is just a comment 
            if line[0] == "#":
                continue 

            comm += line 

            # case: end of command
            if line[-1] == ".": 
                stat = False 

        comm = comm.strip(".")
        self.cmdlines.append(comm)
        return comm 

    # NOTE: does not check for correctness of input 
    def load_cmd_lines(self,cmdlines):
        self.cmdlines.extend(cmdlines)
        return

    def process_command(self): 
        if len(self.cmdlines) == 0: 
            print("NO COMMANDS IN QUEUE")
            return 

        c = self.cmdlines.pop(0)
        c = c.split(" ")
        c = [c_.strip() for c_ in c]
        self.commond = [c_ for c_ in c if len(c_) > 0]
        return self.process_command_()

    def process_command_(self):

        if self.commond[0] == "set": 
            return self.SET_proc(self.commond) 
        return self.MRCOW_proc(self.commond) 

    def SET_proc(self,splitstr_cmd):
        assert splitstr_cmd[0] == "set" 
        n = splitstr_cmd[1] 
        assert splitstr_cmd[2] == "=" 

        self.vartable[n] = self.MRCOW_proc(splitstr_cmd[3:]) 
        return n, self.vartable[n]
    
    """
    method for handling `make`,`run`,`convert`,`open`,`write` processes. 
    """
    def MRCOW_proc(self,splitstr_cmd): 

        if splitstr_cmd[0] == "make":
            return MAKE_proc(splitstr_cmd,self.vartable) 
        elif splitstr_cmd[0] == "run": 
            return RUN_proc(splitstr_cmd,self.vartable)
        elif splitstr_cmd[0] == "convert": 
            return CONVERT_proc(splitstr_cmd,self.vartable)
        elif splitstr_cmd[0] == "open": 
            return OPEN_proc(splitstr_cmd)
        elif splitstr_cmd[0] == "write":
            return WRITE_proc(splitstr_cmd,self.vartable)
        assert False 