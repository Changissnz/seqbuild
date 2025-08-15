import os 

COMM_LANG_GUIDE_FILE = "face/sample_script/commlang_guide.md"

class CLGuideParser:

    def __init__(self): 
        assert os.path.exists(COMM_LANG_GUIDE_FILE) 

        self.file_obj = open(COMM_LANG_GUIDE_FILE,'r') 
        self.file_obj.seek(0,os.SEEK_END) 
        self.file_end = self.file_obj.tell()
        self.file_obj.seek(0)  

        self.structures = dict() # str -> [parameter info] 
        self.keywords = dict() # level -> [keywords] 
        # command -> usage/description -> info. 
        self.command_forms = dict() # command -> usage 

        self.current_line = None 
        return 

    def process(self): 
        self.process_structures() 
        self.process_keywords() 
        self.process_command_forms() 
        return

    def process_structures(self): 
        self.go_to_section("structures") 

        while self.file_obj.tell() != self.file_end:  
            q = self.next_structure()
            if type(q) == type(None):
                break 
            
            n0,n1,n2 = q 
            self.structures[n0] = dict() 
            self.structures[n0]["init"] = n1 
            self.structures[n0]["run"] = n2 

    """
    """
    def next_structure(self): 
        
        def full_parse_line(d,qline_,index): 
            
            stat = True
            if len(qline_) > 0: 
                stat = qline_[-1] != ")"
            
            while stat: 
                line = self.file_obj.readline().strip() 
                qline_ = line 

                if len(qline_) == 0: continue 
                d[index][-1] += qline_
                stat = qline_[-1] != ")" 

        def parse_section(d,s=0):
            c = 0 
            while True: 
                position = self.file_obj.tell() 

                line = self.file_obj.readline().strip()
                if len(line) == 0: 
                    continue

                # case: start of section `instantiation`
                if line == "- instantiation:": 
                    continue 

                # case: start of section `run` 
                if line == "- `run with` parameter:": 
                    # end of `instantiation`
                    if s == 0: 
                        self.file_obj.seek(position) 
                        break 

                    continue 

                # case: end of structure 
                if line[:3] == "[o]": 
                    self.file_obj.seek(position)
                    break 

                # case: end of structure 
                if line[:2] == "##": 
                    self.file_obj.seek(position) 
                    break 

                # case: part of instantiation description
                if line[0] == "*":
                    qline = line.strip("*")
                    qline = qline.strip() 
                    d[c].append(qline)
                    full_parse_line(d,qline,c)
                    c += 1 

                    position = self.file_obj.tell() 
                    line = self.file_obj.readline().strip() 
                    
                    
                    if len(line) == 0: 
                        continue 

                    if line[0] == "*": 
                        c -= 1 
                        self.file_obj.seek(position) 
                # case: start of instantiation description 
                else: 
                    qline = line[2:].strip()
                    d[c] = [qline]
                    full_parse_line(d,qline,c) 
                    
        
        # get the starting line of the next structure 
        line = ""
        stat = False
        while self.file_obj.tell() != self.file_end:  
            position = self.file_obj.tell() 
            line = self.file_obj.readline().strip()

            if line[:3] == "[o]": 
                stat = True
                break 

            if line[:2] == "##": 
                self.file_obj.seek(position) 
                break 

        if not stat: 
            return None 

        line = line[3:].strip() 
        struct_name = line 

        instant_dict = dict()
        run_dict = dict() 
        parse_section(instant_dict,0)
        parse_section(run_dict,1) 

        return struct_name,instant_dict,run_dict 

    def process_keywords(self): 
        self.go_to_section("keywords")

        def process_s12(sx):
            while self.file_obj.tell() != self.file_end:  
                line = self.file_obj.readline() 
                line = line.strip()
                if len(line) == 0: break
                keyword = line[3:].strip() 
                sx |= {keyword} 

        def process_s3(dx): 
            keyword = None 
            while self.file_obj.tell() != self.file_end:  
                line = self.file_obj.readline() 
                line = line.strip()
                if len(line) == 0: break

                if line[:3] == "[-]": 
                    i1 = line.index("`") 
                    sline = line[i1+1:] 
                    i2 = sline.index("`") 
                    keyword = sline[:i2] 
                    dx[keyword] = [] 
                elif line[0] == "*": 
                    kw2 = line[1:].strip() 
                    dx[keyword].append(kw2)
                else: 
                    assert False 

        # process primary,secondary,tertiary in 
        # that order 
        headlines = ["- Primary","- Secondary","- Tertiary"] 

        primary,secondary = set(),set() 
        tertiary = dict()

        self.search_for_line(headlines[0])
        process_s12(primary) 

        self.search_for_line(headlines[1]) 
        process_s12(secondary) 

        self.search_for_line(headlines[2]) 
        process_s3(tertiary) 

        self.keywords[0] = primary 
        self.keywords[1] = secondary
        self.keywords[2] = tertiary
        return primary,secondary,tertiary
    
    def process_command_forms(self): 
        self.go_to_section("command forms") 

        stat = True 

        while stat: 
            c = self.next_command_form()
            stat = not type(c) == type(None) 
            if not stat: continue 

            kw,u,d = c 
            self.command_forms[kw] = (u,d) 
        

        '''
[+] write  
[-] usage  
```
write <object> to <file_object>.  
``` 
[-] description  
writes an object loaded in program memory to a file object. 
        '''
    def next_command_form(self):
        halt_line = "## Interface Layout" 
        line = self.search_for_line("[+]",halt_line=halt_line,is_prefix_search=True)

        if type(line) == type(None): 
            return None 

        keyword = line[3:].strip() 

        # get usage
        usage_ = self.search_for_line("[-] usage",halt_line=halt_line,is_prefix_search=True) 
        usage = "" 
        c = 0 
        while True: 
            if c == 2: 
                break 

            line = self.file_obj.readline().strip() 
            if line == "```": 
                c += 1 
                continue
            usage += line + "\n" 

        # get description
        self.search_for_line("[-] description",halt_line=halt_line,is_prefix_search=True) 
        description = "" 
        while True: 
            line = self.file_obj.readline().strip()
            if len(line) == 0: continue 
            if line[0] == "-": break 
            description += line + "\n" 
        return keyword,usage,description 
            
    def go_to_section(self,section): 
        assert section in {"structures","keywords","command forms"} 
        if section == "command forms": section = "command Forms" 

        section_str = "## " + section[0].upper() + section[1:]
        self.search_for_line(section_str,None,False)

    def search_for_line(self,line_,halt_line=None,is_prefix_search:bool=False):
        while self.file_obj.tell() != self.file_end:  
            line = self.file_obj.readline() 
            line = line.strip()

            if is_prefix_search: 
                if line[:len(line_)] == line_: 
                    return line 

            if line == line_: 
                return line 

            if line == halt_line: 
                return None 

        return None 

    def close(self): 
        self.file_obj.close() 
        del self 