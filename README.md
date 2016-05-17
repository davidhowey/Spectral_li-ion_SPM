About _Spectral li-ion SPM_
===========================

_Spectral li-ion SPM_ is a MATLAB code that solves the so-called 
lithium-ion battery Single Particle Model (SPM) using spectral numerical 
methods.
The SPM is an electrochemical model describing lithium transport, reaction 
kinetics and thermodynamics in lithium-ion batteries. 
_Spectral li-ion SPM_ consists of the SPM coupled to a bulk thermal model 
describing the evolution of battery temperature. 
The SPM is an approximation to the electrochemical pseudo-two dimensional 
lithium-ion battery model where electrolyte transport limitations are 
neglected. 
Therefore the SPM is only valid for relatively low currents, up to 1C or 2C 
depending on the battery design.
The diffusion partial-differential equations of the SPM are discretised in 
space using an efficient spectral numerical method: Chebyshev orthogonal 
collocation. 
This implementation of the SPM is similar to our implementation of the more 
complex pseudo-two dimensional battery model which is discussed in our 
paper:

Bizeray A.M., Zhao S., Duncan S. R., and Howey D. A., 
“Lithium-ion battery thermal-electrochemical model-based state 
estimation using orthogonal collocation and a modified extended Kalman 
filter”, Journal of Power Sources, vol. 296, pp. 400-412, 2015. 
[Publisher copy][6] and [Open access pre-print][7].

If you use _Spectral li-ion SPM_ in your work, please cite our paper.
This code has been developed at the Department of Engineering Science of 
the University of Oxford. 
For information about our lithium-ion battery research, 
visit the [Howey Research Group][2] website.
If you are interested in our energy research, 
check out our research group website [Energy and Power Group][1]. 

For more information and comments, please contact 
[david.howey@eng.ox.ac.uk][5].

Requirements
============
You will need MATLAB to run this code. This code has been developed and 
tested in MATLAB R2015b and should work with later versions. 
Although it has not been tested with earlier MATLAB releases, it should 
also work with no or minor modifications.

You will also need the [_MATLAB Differentiation Matrix Suite_][4] developed 
by Weideman and Reddy to use our code (see "Installation"). The _MATLAB 
Differentiation Matrix Suite_ is available for free on the 
[MathWorks website][3] and is used to compute the Chebyshev orthogonal 
collocation differentiation matrices.
More details can be found on their [website][4] and [paper][10]: 
JAC Weideman, SC Reddy, A MATLAB differentiation matrix suite, 
ACM Transactions of Mathematical Software, 
Vol 26, pp 465-519 (2000).
 
Installation
============
## Step 1 - Installing _Spectral li-ion SPM_ ##
###Option 1 - Downloading a .zip file###
Download a .zip file of the code at:

[https://github.com/adrienBizeray/Spectral_li-ion_SPM/archive/master.zip][8]

Then, unzip the folder in a chosen directory on your computer.

###Option 2 - Cloning the repository with Git###
To clone the repository, you will first need to have [Git][9] installed on 
your computer. Then, navigate to the directory where you want to clone the 
repository in a terminal, and type:
```
git clone https://github.com/adrienBizeray/Spectral_li-ion_SPM.git
```
The folder containing all the files should appear in your chosen directory.

## Step 2 - Installing the _MATLAB Differentiation Matrix Suite_ ##
1. Go to the MathWorks website at 
http://uk.mathworks.com/matlabcentral/fileexchange/29-dmsuite 
and download the `DMSUITE` folder 
(you may need to register and obtain a Mathworks account,
alternatively you can download the `chebdif.m` file only [here][4]).

2. Unzip the `DMSUITE` folder and put it somewhere on your MATLAB path. 
We suggest to drop it in the `source` sub-folder of the 
_Spectral li-ion SPM_ code and run the following command to add this 
sub-folder to the MATLAB path:
```
addpath(genpath('source'));
```

Getting started
===============
The best way to get started with _Spectral li-ion SPM_ is to have a look at
the [EXAMPLE_constant_current_discharge.m](EXAMPLE_constant_current_discharge.m) file. 
This file runs a 1C constant-current discharge simulation for a LCO 
lithium-ion cell, and also gives some explanation about the single particle
model and its MATLAB implementation using Chebyshev orthogonal collocation.

The [EXAMPLE_constant_current_discharge.m](EXAMPLE_constant_current_discharge.m) 
script was written with MATLAB Publish markups, therefore the best way to 
read this file is to generate a HTML file by using the **MATLAB Publish** 
function instead of simply running the file (this may take a few minutes 
the first time as pictures will be automatically generated).
Alternatively, a PDF version of the file is available in the 
_Spectral li-ion SPM_ folder.

License
=======

This open-source MATLAB code is published under the BSD 3-clause License,
please read the `LICENSE.txt` file for more information.

[1]: http://epg.eng.ox.ac.uk/
[2]: http://users.ox.ac.uk/~engs1053/
[3]: http://uk.mathworks.com/matlabcentral/fileexchange/29-dmsuite
[4]: http://dip.sun.ac.za/~weideman/research/differ.html
[5]: mailto:david.howey@eng.ox.ac.uk
[6]: http://www.sciencedirect.com/science/article/pii/S0378775315300677
[7]: http://arxiv.org/abs/1506.08689
[8]: https://github.com/adrienBizeray/SPM/archive/master.zip
[9]: https://git-scm.com/
[10]:http://dl.acm.org/citation.cfm?id=365727