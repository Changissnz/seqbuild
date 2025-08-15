# About: Comm Lang 

## Description

Comm Lang (shorthand for "communication language") is a custom 
language, designed with simplicity in mind, that wraps around 
the code infrastructure of `seqbuild` library. Users of `seqbuild` 
interface can use Comm Lang to write command scripts. Its set of 
keywords is relatively small, much smaller than that of Python 
standard library. And the commands possible to execute via Comm Lang 
are related to `seqbuild` components and capabilities.

## Structures 

Most of the structures are numerical generators. One primary 
structure that is not a numerical generator, but quite important 
in this library, is `MultiMetric`. Below is a list of these structures, 
as well as their instantiation parameters in Comm Lang. 

[o] lcg 
- instantiation:  

    1. (int|float,int|float,int|float,int|float)  
    * (start value,multiple,additive,modulus)  

[o] lcgv2 
- instantiation: 

    1. (int|float,int|float,int|float,int|float,int|float)  
    * (start value,multiple,additive,start modulus,end modulus)  

    2. (int|float,int|float,int|float,int|float,int|float,int,bool)  
    * (start value,multiple,additive,start modulus,end modulus,memory vector size,graph decomposition)  

[o] lcgv3 
- instantiation: 

    1. (int|float,int|float,int|float,int|float,int|float)   
    * (start value,multiple,additive,start modulus,end modulus)     

    2. (int|float,int|float,int|float,int|float,int|float,function,
      int|float,int|float)  
    * (start value,multiple,additive,start modulus,end modulus,
        inputless generator,super-range start,super-range end)   

    3. (int|float,int|float,int|float,int|float,int|float,function,
      int|float,int|float,bool)  
    * (start value,multiple,additive,start modulus,end modulus,
        inputless generator,super-range start,super-range end,exclude zero trinary)  

    4. (int|float,int|float,int|float,int|float,int|float,function,
      int|float,int|float,bool,bool)  
    * (start value,multiple,additive,start modulus,end modulus,
        inputless generator,super-range start,super-range end,
        exclude zero trinary mode, reflective modification mode)  


[o] mdr 
- instantiation:  

    1. (list)  
    * numerical sequence  

[o] mdrv2 
- instantiation:  

    1. (list)  
    * numerical sequence  

    2. (list,bool)  
    * numerical sequence, exclude negative multiples  

[o] mdrgen 
- instantiation:  

    1. (mdr|mdrv2,function,bool,1|2)  
    * (reference `ModuloDecompRepr`,inputless generator,generative type)  

    2. (mdr|mdrv2,function,bool,1|2,bool,bool,bool,bool,bool)  
    * (reference `ModuloDecompRepr`,inputless generator,generative type,
        row-column switch,selector switch 1,selector switch 2, selector switch 3, 
        input seed in output)  

[o] optri  
- instantiation: 

    1. (int,function,1|2,bool,list)  
    * (integer seed,inputless generator,generative type,add noise,base sequence)  

[o] qval  
- instantiation:  

    1. (list,function,function,function,1|2)  
    * (base sequence,inputless generator,inputless generator,
        inputless generator,adjustment type)  
    * (!,index selector,length outputter,range outputter,!)

[o] pid  

- instantiation:  

    1. (function,function,function,1|2)  
    * (inputless generator,inputless generator,inputless generator,
        adjustment type)  
    * (base generator,frequency generator,length generator,range generator,!)  

[o] multimetric  

- instantiation:  

    1. (list)  
    * (base sequence)  

- `run with` parameter:  

    1. (integer)  
    * (positive integer specifying n-gram)  


## Keywords 

- Primary  
    [-] make  
    [-] run  
    [-] set  
    [-] write  
    [-] open  
    [-] convert  

- Secondary  
    [-] with  
    [-] for  
    [-] iter  
    [-] to  

- Tertiary  
    [-] associated with `open`:  
        * file  
        * seq  
        * obj  
    [-] associated with `convert`:  
        * range  
        * ndim  
        * nvec  
        * tvec  

## Command Forms 

[+] make  
[-] usage  
```
make <structure> with <input1,input2,...,inputJ>. 
```
[-] description  
instantiates a structure.  

-----------------------------------------------------------------

[+] run  
[-] usage   
```
run <structure>.  
run <structure> with <input1,...,inputJ>.  
run <structure> for <positive integer> iter.  
```
[-] description  
calls the structure's main output function. 

-----------------------------------------------------------------

[+] set 
[-] usage
```
set <variable name> = <command statement that produces object>. 
```
[-] description  
loads an object into parser map of variables.  

-----------------------------------------------------------------

[+] open  
[-] usage  
```
open file <filepath>.  
open file <filepath> for seq.  
open file <filepath> for obj.   
```
[-] description  
loads a file into a program object. By default, opens file in 
bytes mode to store objects. Specifying `for seq` opens file in 
regular string mode and `for obj` opens file in bytes mode.  

-----------------------------------------------------------------

[+] write  
[-] usage  
```
write <object> to <file_object>.  
``` 
[-] description  
writes an object loaded in program memory to a file object. 

-----------------------------------------------------------------

[+] convert  
[-] usage 
```
convert G to range.  
convert G to ndim with <positive_integer_sequence>.  
convert G to nvec with <positive_integer>.  
convert G to tvec with <positive_integer>.  
```

-----------------------------------------------------------------

## Typical Commands 

## Auxiliary Commands 
