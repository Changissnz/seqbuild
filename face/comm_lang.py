"""
base file for interface language design and corresponding 
interpretive functions to process commands. 
"""

BASE_PRNG = ["lcg","lcgv2","lcgv3","mdr","idforest","optri",\
             "rch"]

LANG_KEYTERMS = ["make","run","with","set","for","iter","write","to"]  

def MAKE_proc(cmd): 
    parsed_cmd = cmd.split(" ") 
    return 