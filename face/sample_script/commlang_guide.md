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
    * (int|float,int|float,int|float,int|float)  
    * (start value,multiple,additive,modulus)  

[o] lcgv2 
- instantiation: 
    * (int|float,int|float,int|float,int|float,int|float)  
    * (start value,multiple,additive,start modulus,end modulus)  

    * (int|float,int|float,int|float,int|float,int|float,int,bool)  
    * (start value,multiple,additive,start modulus,end modulus,memory vector size,graph decomposition)  

[o] lcgv3 
- instantiation: 
    * (int|float,int|float,int|float,int|float,int|float)  
    * (start value,multiple,additive,start modulus,end modulus)  

    * (int|float,int|float,int|float,int|float,int|float,function,
      int|float,int|float)  
    * (start value,multiple,additive,start modulus,end modulus,
        inputless generator,super-range start,super-range end)   

    * (int|float,int|float,int|float,int|float,int|float,function,
      int|float,int|float,bool)  
    * (start value,multiple,additive,start modulus,end modulus,
        inputless generator,super-range start,super-range end,exclude zero trinary)  

    * (int|float,int|float,int|float,int|float,int|float,function,
      int|float,int|float,bool,bool)  
    * (start value,multiple,additive,start modulus,end modulus,
        inputless generator,super-range start,super-range end,
        exclude zero trinary mode, reflective modification mode)  

[o] mdr 
- instantiation:  

    * (list)  
    * numerical sequence  

[o] mdrv2 
- instantiation:  

    * (list)  
    * numerical sequence  

    * (list,bool)  
    * numerical sequence, exclude negative multiples  

[o] mdrgen 
- instantiation:  

    * (mdr|mdrv2,function,bool,1|2)  
    * (reference `ModuloDecompRepr`,inputless generator,generative type)  

    * (mdr|mdrv2,function,bool,1|2,bool,bool,bool,bool,bool)  
    * (reference `ModuloDecompRepr`,inputless generator,generative type,
        row-column switch,selector switch 1,selector switch 2, selector switch 3, 
        input seed in output)  

[o] optri  
- instantiation: 

    * (int,function,1|2,bool,list)  
    * (integer seed,inputless generator,generative type,add noise,base sequence)  

[o] qval  
- instantiation:  

    * (list,function,function,function,1|2)  
    * (base sequence,inputless generator,inputless generator,
        inputless generator,adjustment type)  
    * (!,index selector,length outputter,range outputter,!)

[o] pid  
- instantiation:  

    * (function,function,function,1|2)  
    * (inputless generator,inputless generator,inputless generator,
        adjustment type)  
    * (base generator,frequency generator,length generator,range generator,!)  

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

## Command Forms 

## Typical Commands 

## Auxiliary Commands 
