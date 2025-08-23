import os 

BASE_COMM_LANG_FOLDER = "face/sample_script/"
COMM_LANG_GUIDE_FILE = BASE_COMM_LANG_FOLDER + "commlang_guide.md"

SAMPLE_USE_FILE_MAP = {"multimetric": BASE_COMM_LANG_FOLDER + "commond_one.txt",\
    "lcg":BASE_COMM_LANG_FOLDER + "commond_one.txt",\
    "mdr":BASE_COMM_LANG_FOLDER + "commond_two.txt",\
    "lcgv3":BASE_COMM_LANG_FOLDER + "commond_three.txt",\
    "convert":BASE_COMM_LANG_FOLDER + "commond_four.txt",\
    "qval":BASE_COMM_LANG_FOLDER + "commond_five.txt",\
    "lcgv2":BASE_COMM_LANG_FOLDER + "commond_six.txt",\
    "mdrv2":BASE_COMM_LANG_FOLDER + "commond_six.txt",\
    "mdrgen":BASE_COMM_LANG_FOLDER + "commond_seven.txt",\
    "optri":BASE_COMM_LANG_FOLDER + "commond_eight.txt",\
    "pid":BASE_COMM_LANG_FOLDER + "commond_nine.txt",\
    "open":BASE_COMM_LANG_FOLDER + "commond_ten.txt",\
    "write":BASE_COMM_LANG_FOLDER + "commond_ten.txt",\
    "span":BASE_COMM_LANG_FOLDER + "commond_13.txt",\
    "make":BASE_COMM_LANG_FOLDER + "commond_14.txt",\
    "merge":BASE_COMM_LANG_FOLDER + "commond_15.txt",\
    "load":BASE_COMM_LANG_FOLDER + "commond_16.txt",\
    "rch":BASE_COMM_LANG_FOLDER + "commond_20.txt",\
    "echo":BASE_COMM_LANG_FOLDER + "commond_21.txt"}

def stringize_CLGuideParser_keywords(clgp): 
    term_map = {0:"Primary",1:"Secondary",2:"Tertiary"}
    s = "\t\t[+] Keywords\n\n"  

    for i in range(3): 
        kx = clgp.keywords[i]
        s += "\t{}\n".format(term_map[i])

        if i != 2:
            kx_ = sorted(kx)  
            kx_ = ["* " + kx__ for kx__ in kx_] 
            s += "\n".join(kx_) 
            s += "\n"
        else: 
            kx_ = sorted(kx.keys()) 
            
            while len(kx_) > 0: 
                kx2 = kx_.pop(0) 
                s += "+ " + kx2 + "\n" 

                v2 = kx[kx2] 
                v2 = ["*\t" + v2_ for v2_ in v2]
                s += "\n".join(v2) 
                s += "\n\n"
        
        s += "\n"
    return s 

def stringize_one_CLGuideParser_structure(clgp,struct_name):
    if struct_name not in clgp.structures: 
        return None 

    s_ = "* " + struct_name 
    s_ += "\n\t\tinit:\n"

    v = clgp.structures[struct_name]

    k2 = len(v['init'])
    for i in range(k2):
        s_ += "\n".join(v['init'][i])
        s_ += "\n\n" 

    k2 = len(v['run'])
    if k2 == 0: 
        s_ += "\n"
        s_ += "\n-/-/-/-----------------/-/-/-\n"
        return s_  


    s_ += "\n"
    s_ += "\n\t\trun:\n"
    for i in range(k2): 
        s_ += "\n".join(v['run'][i])
        s_ += "\n"

    s_ += "\n-/-/-/-----------------/-/-/-\n"
    return s_ 

def stringize_CLGuideParser_structures(clgp): 
    assert type(clgp) == CLGuideParser 
    s = "\t\t[+] Structures\n\n" 
    for k in clgp.structures.keys():        
        s_ = stringize_one_CLGuideParser_structure(clgp,k)
        s += s_ 
    return s 

def stringize_one_CLGuideParser_command_form(clgp,command_name): 
    if command_name not in clgp.command_forms: 
        return None
    
    v = clgp.command_forms[command_name] 
    usage,description = v[0],v[1] 

    s = "* " + command_name + "\n" 
    s += "\t- usage\n" + usage + "\n"
    s += "\t- description\n" + description + "\n"
    s += "\n-/-/-/-----------------/-/-/-\n"
    s += "\n"
    return s 
    
def stringize_CLGuideParser_command_forms(clgp): 

    assert type(clgp) == CLGuideParser 
    s = "\t\t[+] Command Forms\n\n" 

    for k in clgp.command_forms.keys():
        s += stringize_one_CLGuideParser_command_form(clgp,k)

    return s 

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

        self.has_processed = False 
        return 

    def about(self,term):
        ktype_dict = {0:"primary ", 1:"secondary ",2:"tertiary ",-1:"not a "} 
        ktype = self.keyword_type(term) 
        s0 = ktype_dict[ktype] + " keyword\n" 

        s = "" 
        if term in self.structures: 
            s0 += "\t* structure\n"
            s = stringize_one_CLGuideParser_structure(self,term) 
        elif term in self.command_forms:
            s0 += "\t* command\n"
            s = stringize_one_CLGuideParser_command_form(self,term)
            
        if type(s) == type(None): 
            s = "" 
        s = s0 + "\n-/-/-/-----------------/-/-/-\n" + s

        if term in SAMPLE_USE_FILE_MAP: 
            fi_obj2 = open(SAMPLE_USE_FILE_MAP[term],"r") 
            lines = fi_obj2.readlines() 
            s = s + "\n\t\tsample use\n---------------\n"
            s += ''.join(lines)
            fi_obj2.close() 
        
        return s 
 
    def keyword_type(self,term):
        if term in self.keywords[0]:
            return 0 
        
        if term in self.keywords[1]: 
            return 1 

        q = self.keywords[2] 
        x = [] 
        for v_ in q.values(): 
            x.extend(list(v_)) 
        if term in x: return 2  

        return -1 

    def __str__(self):
        assert self.has_processed

        s1 = stringize_CLGuideParser_keywords(self) 
        s2 = stringize_CLGuideParser_structures(self) 
        s3 = stringize_CLGuideParser_command_forms(self) 
        return s1 + "_" * 75 + "\n\n" \
            + s2 + "_" * 75 + "\n\n" + s3 

    """
    main method 
    """
    def process(self): 
        self.process_structures() 
        self.process_keywords() 
        self.process_command_forms()
        self.has_processed = True  
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
        self.search_for_line(section_str,None,True) 

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