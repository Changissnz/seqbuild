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

    def go_to_section(self,section): 
        assert section in {"structures","keywords","command forms"} 
        section_str = "## " + section[0].upper() + section[1:]

        while self.file_obj.tell() != self.file_end:  
            line = self.file_obj.readline() 
            line = line.strip()

            if line == section_str: 
                break

    def close(self): 
        self.file_obj.close() 
        del self 