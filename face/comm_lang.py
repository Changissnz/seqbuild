"""
base file for interface language design and corresponding 
interpretive functions to process commands. 
"""
from .cmd_proc import * 

class CommLangParser:

    def __init__(self,filepath:str):
        self.fp = filepath
        self.file_obj,self.file_end = None,None
        self.line_index = 0 
        self.preproc_file() 

        self.cmdlines = []
        self.show_commands = [] 
        self.help_commands = [] 
        self.extra_commands = [] 

        self.cmd_errors = []

        self.commond = None 
        self.vartable = dict() 

    def __str__(self):

        def map_str(s): 
            if s != "method": return s 
            return "generator" 

        dx = ""
        for k,v in self.vartable.items(): 
            s = str(k) 
            s = s + "\t" + map_str(parse_object_to_str(v))
            dx += s + "\n" 
        return dx 

    def preproc_file(self): 
        if type(self.fp) == type(None): 
            return 

        self.file_obj = open(self.fp,"r") 
        self.file_obj.seek(0,os.SEEK_END) 
        self.file_end = self.file_obj.tell()
        self.file_obj.seek(0)  
        self.line_index = 0 
        self.cmd_errors = [] 

    def reload_file(self,filepath:str): 
        self.fp = filepath 
        self.preproc_file() 

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

    """
    used to process all remaining commands from `cmdline` cache. 
    """
    def process_cmdlines(self): 
        while len(self.cmdlines) > 0: 
            self.process_command()

    def load_next_command(self):
        comm = "" 
        stat = True 

        while stat: 
            if self.file_obj.tell() != self.file_end:
                line = self.file_obj.readline()
                self.line_index += 1 
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

        c_ = self.cmdlines.pop(0)
        c = c_.split(" ")
        c[-1] = c[-1].rstrip(".")

        # case: `show` command, do nothing 
        if c[0] == "show": 
            self.show_commands.append(c_) 
            return 

        # case: `help` command, do nothing 
        if c[0] == "help":
            self.help_commands.append(c_) 
            return 

        if c[0] == "load": 
            qs = LOAD_proc(c) 
            rx = deepcopy(self.cmdlines)
            self.cmdlines.clear() 
            while len(qs) > 0:
                qs_ = qs.pop(-1).strip()
                self.cmdlines.insert(0,qs_) 
            self.process_cmdlines()
            self.cmdlines = rx 
            return 

        if c[0] == "qualtest":
            self.extra_commands.append(c)
            return 

        if c[0] == "chaintest":
            self.extra_commands.append(c) 
            return 

        self.commond = [c_ for c_ in c if len(c_) > 0]
        return self.process_command_()

    def process_command_(self):
        if len(self.commond) == 0: 
            return 
            
        if self.commond[0] == "set": 
            try: 
                q = self.SET_proc(self.commond) 
                return q 
            except: 
                self.log_error(self.commond) 
                return 

        return self.MRCOW_proc(self.commond) 

    def SET_proc(self,splitstr_cmd):
        assert splitstr_cmd[0] == "set" 
        n = splitstr_cmd[1] 
        assert n not in LANG_KEYTERMS 
        assert splitstr_cmd[2] in {"=","span"}  

        # case: primary use case of SET 
        if splitstr_cmd[2] == "=": 
            self.vartable[n] = self.MRCOW_proc(splitstr_cmd[3:]) 
            return n, self.vartable[n] 
        
        # case: secondary use case of SET 
        assert n in self.vartable 
        gen = MAIN_method_for_object(self.vartable[n])

        assert splitstr_cmd[3] == "to"  

        rx = splitstr_cmd[4].strip(".").split(",")
        assert len(rx) == 2

        f0 = float(rx[0]) if ("." in rx[0] \
            or "E" in rx[0]) else int(rx[0]) 
        f1 = float(rx[1]) if ("." in rx[1] \
            or "E" in rx[1]) else int(rx[1])

        if type(f0) == float or type(f1) == float: 
            f0,f1 = float(f0),float(f1) 

        try:
            x = wrap_ranged_modulo_over_generator(gen,(f0,f1))  
            self.vartable[n] = x 
        except: 
            print("cannot set span for non-generator.") 
            self.log_error(splitstr_cmd) 
        return 

    
    """
    method for handling `make`,`run`,`convert`,`open`,`write` processes. 
    """
    def MRCOW_proc(self,splitstr_cmd): 
        try: 
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
            elif splitstr_cmd[0] == "merge": 
                return MERGE_proc(splitstr_cmd,self.vartable) 
            elif splitstr_cmd[0] == "read": 
                return READ_proc(splitstr_cmd)
        except: 
            pass 

        self.log_error(splitstr_cmd) 
        return None 

    def log_error(self,splitstr_cmd):
        q = ' '.join(splitstr_cmd)
        s = '[file line:' + str(self.line_index) + ']' + q 
        self.cmd_errors.append(s) 