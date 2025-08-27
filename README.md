# About: seqbuild 

Sequentete...sequentete, my son. 

*Work-in-progress* 

Project is being developed in conjunction with 
that of dependency library, `morebs2`. Started out 
as a research project, and I am starting to wrap 
ends to it before the project branches off into 
development hell terrain. 

And here is this, because and for the 
kleptos, and because not specifying 
copyright encourages intellectual 
property theft (I may assume IP laws 
are still honored). Remember that all 
software is proprietary unless otherwise 
specified by a FOSS license. 

*Sorry that I had to put this poison pill provision into this project.*
*That said, the free market has never been that free with research and*
*emerging technologies.*

Copyright 2025 Richard Pham 

### Brief Description (Current State)  

Primary use is research of pseudo-random number generators (PRNG). 
Comm Lang is the application programming interface that `seqbuild` uses, 
and users can express commands to gauge and generate numerical sequences 
in this language. When user runs `main.py` (the main program), they are 
accessing a semi-autonomous PRNG system for their personal use. They are 
responsible for: 
- initializing generators,  
- gauging generators and numerical sequences,  
- tweaking generators and numerical sequences.  

Autonomous capabilities for `seqbuild` to perpetually output values classified 
as "random", according to certain metrics, have not been implemented. The capabilities 
of reading generators from and writing generators to files have also not been 
implemented. This is due to the essential role of the `pickle` library. The 
`pickle` library does not support the storage of Python objects that are closures 
or contain closures. A closure in the Python programming language is a function nested 
inside another function; there are variants with "lambda","method", and "function". 
The `json` library also does not store objects of or related to closures. Without a 
convenient methodology for reading generators from and writing generators to files, 
it is easy to see why these capabilities have not been implemented, however important 
they are.  

### Update: 8/27/25  

I have decided to add a new folder, `performance_marks`, that contain files named  
`fYxx.txt`, along with their corresponding randomness test results `perfY.txt`. The  
tests are done by an implementation of MNIST Randomness Testing Suite. And this implementation 
can be found at the link below.  

`https://github.com/stevenang/randomness_testsuite.git`. 

These test results demonstrate the challenge of ?defining the concept? of random. 

### Update: 8/26/25  

The main program now includes an `encrypt` Comm Lang command. In cases where `np.nan`, 
`np.inf`, or `None` is the ciphertext element, there is not an easy way to map to the 
plaintext character from it. This particular `encrypt` command can be used for 
one-way encryption.

### Update: 8/22/25  

There are still unfinished features in this project. Development will no 
longer proceed at a regular frequency. Rolling developments, of indefinite 
timing, will be pushed to the project.

### Update: 8/14/25

Development is halfway done. The `seqbuild` project now has an interface. 
Here are the steps to using it. 
- go into a terminal prompt and enter in these commands: 
```
> git clone https://www.github.com/changissnz/seqbuild.git 
> cd seqbuild 
> pip install -r requirements.txt
> python3 main.py 
```

Here is my developer note that pops up every time the `main.py` file 
is run. 

![Local Image](2025-08-14__developer_note.png)  

I have planned for 1 full week of additional development, starting today. 
This will be AGILE end sprint for beta version of this product. 