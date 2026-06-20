from .comm_lang import * 
from morebs2.numerical_generator import sign_preserving_modulo
import time 

DEFAULT_TBCLF_TIMESTAMP_SEEDSIZE_RANGE = [2,7] 

DEFAULT_TBCLF_LCG_PARAMETER_RANGE = [-10** 4 - 616/617, 10 ** 5 + 666/667]
DEFAULT_TBCLF_GEN_LCGBUNCH_SIZE_RANGE = [2,7] 
DEFAULT_TBCLF_GEN_ZERO_DEFAULT = 59163 + 56/57 
DEFAULT_TBCLF_GEN_VECTOR_LENGTH_RANGE = [7,19] 

DEFAULT_TBCLF_MODULODECOMP_MAX_ABSMULT_RANGE = [5,50] 

DEFAULT_TBCLF_GEN_RCH_QCAP_RANGE = [321,987]
DEFAULT_TBCLF_GEN_RCH_NODESIZE_RANGE = [6,12]

DEFAULT_TBCLF_GEN_SSINO_NODESIZE_RANGE = [10,42] 

DEFAULT_TBCLF_GEN_IOMAPS_SIZE_RANGE = [2,7] 

DEFAULT_TBCLF_GEN_IDF_CACHE_SIZE_RANGE = [100,455] 

# NOTE: generator may use vectors of length not in DEFAULT_TBCLF_GEN_VECTOR_LENGTH_RANGE, depending 
#       on the parameters `vector_files` and `comm_lang_files`. 
class TimeBasedCommLangFileGenerator: 

    def __init__(self,filepath:str,base_gen_name:str,use_prng_for_prng_pr:float,vector_files=[],\
        comm_lang_files=[]):

        assert type(filepath) == str == type(base_gen_name) 
        assert " " not in filepath
        assert is_alphanumeric(base_gen_name)
        assert base_gen_name != "DEFAULT" 
        assert type(use_prng_for_prng_pr) in {int,float}  
        assert 0. <= use_prng_for_prng_pr <= 1.0 
        assert type(vector_files) == type(comm_lang_files) == list 

        self.filepath = filepath 
        self.base_gen_name = base_gen_name 
        self.gen_count = 0 

        self.use_prng_for_prng_pr = use_prng_for_prng_pr
        self.vfiles = vector_files
        self.vf_index = 0 
        self.vf_index_2 = 0 
        self.cl_files = comm_lang_files 
        
        self.t = time.time()

        q = modulo_in_range(int(self.t),DEFAULT_TBCLF_TIMESTAMP_SEEDSIZE_RANGE)

        x = -1 
        for _ in range(q -1): 
            self.t += time.time() * x 
            x *= -1 

        self.pdd = PRNGDecimalDelta(self.t)

        # command lines for Comm Lang file to be generated 
        self.cl_lines = [] 
        self.valid_fp = True 
        self.clp = None 
        self.erroneous_files = [] 

        # names of available vectors 
        self.available_vectors = [] 

        try: 
            with open(self.filepath,"w") as f: pass 
        except: 
            self.valid_fp = False 
        
        if not self.valid_fp: return 

        self.clp = CommLangParser(self.filepath)

        self.preproc()
        return

    #---------------------------------------------------------- preprocessing methods 

    def errors_exist(self): 
        if len(self.erroneous_files) > 0: 
            return True 

        if len(self.clp.cmd_errors) > 0: 
            return True 
        
    def preproc(self): 
        # check all the vector files 
        for vf in self.vfiles: 
            stat = os.path.isfile(vf) 
            if not stat: 
                self.erroneous_files.append(vf) 
        
        if self.errors_exist(): 
            return 

        # load up every vector into file 
        for vf in self.vfiles: 
            self.load_vector_file(vf) 
            if self.errors_exist(): 
                return 

        # load up every Comm Lang file in order 
        for cf in self.cl_files: 
            s = "load {}.".format(cf) 
            self.update_CL_file([s])

            # invalid command file 
            if len(self.clp.cmd_errors) > 0: 
                self.erroneous_files.append(cf) 
                break 

        S = self.clp.object_list("sequence",True)
        self.available_vectors.extend([s[0] for s in S]) 
        return 

    def load_vector_file(self,vfile): 
        vname = self.next_vector_name() 

        G = self.fetch_prg(True) 
        vlength = modulo_in_range(int(G()),DEFAULT_TBCLF_GEN_VECTOR_LENGTH_RANGE)

        s = ["set {} = read {} for {} iter.".format(vname,vfile,vlength)] 
        self.update_CL_file(s) 
        return vname

    def next_generator_name(self): 
        s = self.base_gen_name + "_" + str(self.gen_count)
        self.gen_count += 1 
        return s 

    def next_vector_name(self): 
        s = self.base_gen_name + "_" + "vec_" + str(self.vf_index) 
        self.vf_index += 1 
        return s 

    #----------------------------------- generate LCG* commands 

    def generate_CL_LCG_bunch(self):

        bunch_size = modulo_in_range(int(next(self.pdd)),DEFAULT_TBCLF_GEN_LCGBUNCH_SIZE_RANGE)

        return -1 

    def generate_CL_LCG(self): 

        prg = self.fetch_prg(True) 
        r = self.output_n_values(n=4,prg=prg,retry_nonzero=3,nonzero_indices=set(3))

        gen_name = self.next_generator_name() 
        s = "set {} = make lcg with {},{},{},{}.".format(\
            self.base_gen_name,r[0],r[1],r[2],r[3]) 

        return [s],gen_name 

    def generate_CL_LCGv2(self): 

        prg = self.fetch_prg(True)
        r = self.generate_base5_LCG_numbers(prg)            

        gen_name = self.next_generator_name() 
        s = "set {} = make lcgv2 with {},{},{},{},{}.".format(\
            gen_name,r[0],r[1],r[2],r[3],r[4]) 

        return [s],gen_name

    def generate_CL_LCGv3(self): 
        prg_ = self.fetch_prg(True)
        base_numbers = self.generate_base5_LCG_numbers(prg_)

        prg = self.fetch_accessory_prg_varname() 

        gen_type = round(next(self.pdd)) % 5 

        gen_name = self.next_generator_name() 
        if gen_type == 0: 
            s = "set {} = make lcgv3 with {},{},{},{},{},{}.".format(\
                gen_name,base_numbers[0],base_numbers[1],base_numbers[2],\
                base_numbers[3],base_numbers[4],prg) 
            return [s],gen_name

        b = int(next(self.pdd)) % 2 
        if gen_type == 1: 
            s = "set {} = make lcgv3 with {},{},{},{},{},{},{}.".format(\
                gen_name,base_numbers[0],base_numbers[1],base_numbers[2],\
                base_numbers[3],base_numbers[4],prg,b) 
            return [s],gen_name

        sr_start,sr_end = sorted([prg_(),prg_()]) 
        sr_start,sr_end = self.output_n_values(n=2,prg=prg_,\
            retry_nonzero=3,nonzero_indices={1}) 
        if sr_start == sr_end: 
            sr_end += DEFAULT_TBCLF_GEN_ZERO_DEFAULT

        Q = sorted([base_numbers[3],base_numbers[4],sr_start,sr_end]) 
        base_numbers[3],base_numbers[4] = Q[1],Q[2] 
        sr_start,sr_end = Q[0],Q[3]

        if gen_type == 2: 
            s = "set {} = make lcgv3 with {},{},{},{},{},{},{}.".format(\
                gen_name,base_numbers[0],base_numbers[1],base_numbers[2],\
                base_numbers[3],base_numbers[4],prg,sr_start,sr_end,b)  
            return [s],gen_name 

        ts_range_start,ts_range_end = self.generate_int_range(prg_,DEFAULT_LCGV3_TERNARY_SIZE_RANGE)
        if gen_type == 3:  

            s = "set {} = make lcgv3 with {},{},{},{},{},{},{},{},{}.".format(\
                gen_name,base_numbers[0],base_numbers[1],base_numbers[2],\
                base_numbers[3],base_numbers[4],prg,sr_start,sr_end,ts_range_start,ts_range_end,b) 
            return [s],gen_name 

        tdt_start,tdt_end = self.generate_int_range(prg_,DEFAULT_LCGV3_TERNARY_DELTA_TIMESTAMP_RANGE)
        s = "set {} = make lcgv3 with {},{},{},{},{},{},{},{},{}.".format(\
            gen_name,base_numbers[0],base_numbers[1],base_numbers[2],\
            base_numbers[3],base_numbers[4],prg,sr_start,sr_end,ts_range_start,\
            ts_range_end,tdt_start,tdt_end,b) 

        return [s],gen_name 

    #------------------------------------------- generate MDR commands 

    def generate_CL_mdr(self,mdr_type):

        assert mdr_type in {1,2} 

        # retrieve a vector 
        vname = self.fetch_vector_varname() 

        # get max absolute fitting multiple 
        G = self.fetch_prg(True,True)
        max_absmult = modulo_in_range(int(G()),DEFAULT_TBCLF_MODULODECOMP_MAX_ABSMULT_RANGE) 

        gen_name = self.next_generator_name() 
        
        if mdr_type == 1: 
            s = "set {} = make mdr with {},{}.".format(gen_name,vname,max_absmult)
        else: 
            b = int(G()) % 2 
            s = "set {} = make mdrv2 with {},{},{}.".format(gen_name,vname,b,max_absmult)

        return [s],gen_name 

    def generate_CL_mdrgen(self):  
        
        # generate the `ModuloDecompRepr`
        mdr_type = modulo_in_range(int(next(self.pdd)),[1,3]) 
        s0,g = self.generate_CL_mdr(mdr_type)
        prg = self.fetch_accessory_prg_varname() 

        exclude_neg = int(next(self.pdd)) % 2 
        generative_type = modulo_in_range(int(next(self.pdd)),[1,3])

        q = int(next(self.pdd)) % 2  

        gen_name = self.next_generator_name()

        # declaration #1 
        if q: 
            s = "set {} = make mdrgen with {},{},{},{}.".format(gen_name,g,prg,exclude_neg,generative_type)

        else: 
            G = self.fetch_prg(True,True) 
            G = prg__single_to_int(G) 
            G = wrap_ranged_modulo_over_generator(G,[0,2]) 
            B = prg__single_to_nvec(G,5)()
            B = [int(b) for b in B]

            s = "set {} = make mdrgen with {},{},{},{},{},{},{},{},{}.".format(\
                gen_name,g,prg,exclude_neg,generative_type,B[0],B[1],B[2],B[3],B[4]) 

        return s0 + [s],gen_name 

    #-------------------------------------------- generate commands for `optri`,`qval`

    def generate_CL_optri(self): 
        
        m = int(modulo_in_range(next(self.pdd),DEFAULT_TBCLF_LCG_PARAMETER_RANGE)) 

        prg = self.fetch_accessory_prg_varname() 

        generative_type = modulo_in_range(int(next(self.pdd)),[1,3])
        add_noise = int(next(self.pdd)) % 2 

        base_sequence = self.fetch_vector_varname() 

        gen_name = self.next_generator_name() 
        s = "set {} = make optri with {},{},{},{},{}.".format(\
            gen_name,m,prg,generative_type,add_noise,base_sequence) 

        return [s],gen_name 

    def generate_CL_qval(self): 
        base_sequence = self.fetch_vector_varname() 
        V = self.clp.vartable[base_sequence]
        l = len(V) 

        prg = self.fetch_accessory_prg_varname() 
        #ndim_gen_name = self.next_generator_name() 
        #s0 = "set {} = convert {} to ndim with {},{}.".format(ndim_gen_name,prg,l,l) 
        s0,ndim_gen_name = self.convert_prg_output(prg,\
            "ndim with {},{}".format(l,l)) 

        lout_gen_name_ = self.fetch_accessory_prg_varname()
        #lout_gen_name = self.next_generator_name() 
        #s1 = "set {} = convert {} to int.".format(lout_gen_name,lout_gen_name_) 
        s1,lout_gen_name = self.convert_prg_output(lout_gen_name_,\
            "int") 

        prg1 = self.fetch_accessory_prg_varname() 
        #rout_gen_name = self.next_generator_name() 
        #s2 = "set {} = convert {} to range.".format(rout_gen_name,prg1) 
        s2,rout_gen_name = self.convert_prg_output(prg1,\
            "range") 

        mode = modulo_in_range(int(next(self.pdd)), [1,3])

        qval_gen_name = self.next_generator_name() 
        s3 = "set {} = make qval with {},{},{},{},{}.".format(\
            qval_gen_name,base_sequence,ndim_gen_name,lout_gen_name,rout_gen_name,mode) 

        return [s0,s1,s2,s3], qval_gen_name

    #----------------------------------------------- generate commands for `rch`,`pid`

    def generate_CL_rch(self): 

        q = int(next(self.pdd)) % 3 
        num_nodes = modulo_in_range(int(next(self.pdd)),DEFAULT_TBCLF_GEN_RCH_NODESIZE_RANGE)

        prg = self.fetch_accessory_prg_varname() 
        G = self.fetch_prg(True,True) 
        mutation_rate = round(prg_decimal(G,[0.2,0.8]),5) 
        gen_name = self.next_generator_name() 

        if q == 0: 
            s = "set {} = make rch with {},{},{}.".format(gen_name,num_nodes,prg,mutation_rate)
            return [s],gen_name

        queue_capacity = modulo_in_range(int(G()),DEFAULT_TBCLF_GEN_RCH_QCAP_RANGE)

        if q == 1: 
            s = "set {} = make rch with {},{},{},{}.".format(gen_name,num_nodes,prg,mutation_rate,queue_capacity)
            return [s],gen_name

        r0,r1 = self.output_n_values(2,G,retry_nonzero=3) 
        r0 = sign_preserving_modulo(r0,DEFAULT_TBCLF_GEN_ZERO_DEFAULT)
        r1 = sign_preserving_modulo(r1,DEFAULT_TBCLF_GEN_ZERO_DEFAULT)

        if r0 == r1: 
            r1 += DEFAULT_TBCLF_GEN_ZERO_DEFAULT / 4 

        r0,r1 = sorted([r0,r1]) 

        s = "set {} = make rch with {},{},{},{},{},{}.".format(gen_name,num_nodes,prg,mutation_rate,queue_capacity,\
            r0,r1) 
        return [s],gen_name 

    def generate_CL_pid(self): 
        prg1 = self.fetch_accessory_prg_varname() 
        prg2 = self.fetch_accessory_prg_varname() 
        prg3 = self.fetch_accessory_prg_varname() 
        prg4 = self.fetch_accessory_prg_varname()

        #name1 = self.next_generator_name() 
        #s1 = "set {} = convert {} to nvec with 2.".format(name1,prg1)
        s1,name1 = self.convert_prg_output(prg1,"nvec with 2") 

        #name2 = self.next_generator_name() 
        #s2 = "set {} = convert {} to int.".format(name2,prg2) 
        s2,name2 = self.convert_prg_output(prg2,"int") 

        #name3 = self.next_generator_name() 
        #s3 = "set {} = convert {} to int.".format(name3,prg3) 
        s3,name3 = self.convert_prg_output(prg3,"int") 


        #name4 = self.next_generator_name() 
        #s4 = "set {} = convert {} to range.".format(name4,prg4) 
        s4,name4 = self.convert_prg_output(prg4,"range") 


        G = self.fetch_prg(True,True)  
        adjustment_type = modulo_in_range(int(G()),[1,3])

        gen_name = self.next_generator_name() 
        s = "set {} = make pid with {},{},{},{},{}.".format(gen_name,\
            name1,name2,name3,name4,adjustment_type) 

        return [s1,s2,s3,s4,s],gen_name 

    #--------------------------------------------------- generate commands for `ssino`,`idforest`

    def generate_CL_ssino(self): 

        g = self.fetch_prg(True,True) 
        gen_type = int(g()) % 5 

        num_nodes = modulo_in_range(int(g()),DEFAULT_TBCLF_GEN_SSINO_NODESIZE_RANGE) 

        gen1 = self.fetch_accessory_prg_varname() 
        gen2 = self.fetch_accessory_prg_varname() 

        gen_name = self.next_generator_name()
        if gen_type == 0: 
            s = "set {} = make ssino with {},{},{}.".format(gen_name,\
                num_nodes,gen1,gen2) 
            return [s],gen_name 
        
        F = prg__single_to_int(g)
        F = wrap_ranged_modulo_over_generator(F,[0,2]) 
        B = prg__single_to_nvec(F,4)()
        B = [int(b) for b in B]
        ##print("BBB: ",B)
        ##for _ in range(10): print(F())
        if gen_type == 1:
            s = "set {} = make ssino with {},{},{},{}.".format(gen_name,\
                num_nodes,gen1,gen2,B[0]) 
            return [s],gen_name 

        if gen_type == 2:
            s = "set {} = make ssino with {},{},{},{},{}.".format(gen_name,\
                num_nodes,gen1,gen2,B[0],B[1]) 
            return [s],gen_name 

        if gen_type == 3:
            s = "set {} = make ssino with {},{},{},{}.".format(gen_name,\
                num_nodes,gen1,gen2,B[0],B[1],B[2]) 
            return [s],gen_name 

        s = "set {} = make ssino with {},{},{},{}.".format(gen_name,\
            num_nodes,gen1,gen2,B[0],B[1],B[2],B[3]) 
        return [s],gen_name 

    # NOTE: does not cover all possible instantiation routes for `idforest`. 
    def generate_CL_idforest(self): 
        V = self.fetch_vector_varname() 
        s0, ioname = self.generate_CL_iomaps() 
        prg = self.fetch_accessory_prg_varname()

        gen_name = self.next_generator_name() 

        G = self.fetch_prg(True,True)
        gen_type = int(G()) % 4 

        gen_name_1 = self.next_generator_name() 
        s1 = "set {} = convert {} to int.".format(gen_name_1,prg) 

        if gen_type == 0:         
            s = "set {} = make idforest with {},{},{}.".format(\
                gen_name,V,ioname,gen_name_1)  
            return [s0,s1,s],gen_name   
            

        if gen_type == 1: 
            gen_name_2 = self.next_generator_name()
            prg2 = self.fetch_accessory_prg_varname() 
            s2 = "set {} = convert {} to int.".format(gen_name_2,prg2) 

            s = "set {} = make idforest with {},{},{},{}.".format(\
                gen_name,V,ioname,gen_name_1,gen_name_2)  
            return [s0,s1,s2,s],gen_name   
        
        cache_size = modulo_in_range(int(G()),DEFAULT_TBCLF_GEN_IDF_CACHE_SIZE_RANGE)

        if gen_type == 2:         
            s = "set {} = make idforest with {},{},{},{}.".format(gen_name,V,ioname,\
                cache_size,gen_name_1)
            return [s0,s1,s],gen_name 
         

        gen_name_2 = self.next_generator_name()
        prg2 = self.fetch_accessory_prg_varname() 
        s2 = "set {} = convert {} to int.".format(gen_name_2,prg2) 

        s = "set {} = make idforest with {},{},{},{},{}.".format(gen_name,V,ioname,\
            cache_size,gen_name_1,gen_name_2)
        return [s0,s1,s2,s],gen_name 

    #--------------------------------------------------------- choose active PRNG 

    def convert_prg_output(self,prg,convert_str): 
        convert_type = convert_str.split(" ")[0]  
        assert convert_type in GENFORM_CONVERT_TYPES

        #prg = self.fetch_accessory_prg_varname()
        gen_name = self.next_generator_name() 
        s1 = "set {} = convert {} to {}.".format(gen_name,prg,convert_str)  
        return s1,gen_name 

    def fetch_accessory_prg_varname(self): 
        prg = self.fetch_prg(False,False) 
        if type(prg) == type(None): 
            prg = self.default_PRNG__CL_command()
        return prg 

    def fetch_prg(self,get_object:bool,resort_to_default:bool=True): 

        use_prng_for_prng = prg_decimal(self.pdd.__next__,[0.,1.]) < self.use_prng_for_prng_pr
        return self.fetch_prg_(use_prng_for_prng,get_object,resort_to_default)  

    def fetch_prg_(self,use_prng_for_prng,get_object:bool,resort_to_default:bool=True): 
        if not use_prng_for_prng: 
            if not resort_to_default: return None 

            if get_object: 
                return self.pdd.__next__
            return "DEFAULT"

        full = not get_object
        g_list = self.clp.single_output_generator_list() 

        if len(g_list) == 0:
            if resort_to_default:  
                return self.pdd.__next__ if get_object else "DEFAULT"  
            return None 

        t = round(next(self.pdd)) % len(g_list) 
        x = g_list[t]
        if get_object: return MAIN_method_for_object(self.clp.vartable[x])
        return x 

    def fetch_vector_varname(self): 
        # case: no more vectors 
        if len(self.available_vectors) == 0: 
            v = self.load_generated_vector__CL_command()
        # case: next vector in queue 
        else: 
            v = self.available_vectors.pop(0) 

        return v 

    #-------------------------------------- generating Comm Lang commands for special cases 

    """
    used in cases where a PRNG is needed, but not are available. 
    """
    def default_PRNG__CL_command(self): 
        s,gen_name = self.generate_CL_LCGv2() 
        self.update_CL_file(s)
        return gen_name 

    def load_generated_vector__CL_command(self):
        gname = self.fetch_prg(False,False)  
        if type(gname) == type(None): 
            gname = self.default_PRNG__CL_command() 

        vlength = modulo_in_range(int(next(self.pdd)),DEFAULT_TBCLF_GEN_VECTOR_LENGTH_RANGE)
        vname = self.next_vector_name()
        s = "set {} = run {} for {} iter.".format(vname,gname,vlength)
        self.update_CL_file([s]) 
        self.available_vectors.append(vname) 
        return vname 

    def PRNG_seq_to_PRNG__CL_command(self): 
        return -1 

    def generate_CL_iomaps(self):

        g = self.fetch_prg(True,True) 
        num_generators = modulo_in_range(int(g()),DEFAULT_TBCLF_GEN_IOMAPS_SIZE_RANGE) 

        S = [self.fetch_accessory_prg_varname() for _ in range(num_generators)]
        S_ = ",".join(S) 

        gen_name = self.next_generator_name() 
        s = "set {} = make iomaps with {}.".format(gen_name,S_)
        return s,gen_name

    #----------------------------------------------------------- loading commands to Comm Lang cache and writing out to file 

    def update_CL_file(self,command_list):

        self.clp.load_cmd_lines(command_list)  
        self.clp.process_cmdlines() 
        self.cl_lines.extend(command_list)
        return 

    def write_out_to_CL_file(self):
        s = "\n".join(self.cl_lines) 

        with open(self.filepath,"w") as f: 
            f.write(s)
            f.flush()  
        return

    def close(self): 
        self.clp.close() 
        del self 

    #----------------------------------------------------------- auxiliary commands 

    def output_n_values(self,n,prg,retry_nonzero:int,nonzero_indices=set()):
        assert n >= 0 and type(n) == int 
        assert type(nonzero_indices) == set 
        if len(nonzero_indices) > 0: 
            assert 0 <= min(nonzero_indices) <= max(nonzero_indices) < n

        def go_with_default_nonzero(): 
            d = prg_decimal(prg,[0.,1.]) 
            q = 1 
            if d < 0.5: 
                q *= -1 
            return q * DEFAULT_TBCLF_GEN_ZERO_DEFAULT 

        def retry_nonzero_output():  
            for _ in range(retry_nonzero): 
                x = modulo_in_range(prg(),DEFAULT_TBCLF_LCG_PARAMETER_RANGE) 
                x = round(x,5) 
                if x != 0: 
                    return x 
            return go_with_default_nonzero() 

        r = [round(modulo_in_range(prg(),DEFAULT_TBCLF_LCG_PARAMETER_RANGE),5) \
            for _ in range(n)]

        nzi = sorted(nonzero_indices) 

        for n in nzi: 
            if r[n] == 0: 
                r[n] = retry_nonzero_output() 
        return r 

    def generate_base5_LCG_numbers(self,prg): 
        r = self.output_n_values(n=5,prg=prg,\
            retry_nonzero=3,nonzero_indices={4})
        if r[3] == r[4]: 
            r[4] += DEFAULT_TBCLF_GEN_ZERO_DEFAULT

        if r[3] > r[4]: 
            r[3],r[4] = r[4],r[3] 
        return r 

    def generate_int_range(self,prg,R): 
        assert is_valid_range(R,True,False) and R[1] - R[0] > 1 
        S = modulo_in_range(int(prg()),R) 
        E = modulo_in_range(int(prg()),R) 

        if S == E: 
            E = modulo_in_range(E + 1,R)

        S,E = sorted([S,E])
        return S,E 

        
################################################