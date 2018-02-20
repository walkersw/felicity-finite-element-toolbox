How to Install FELICITY
=======================

Installing FELICITY mostly involves unzipping it to a directory, configuring the MATLAB path, and running some tests to verify that it works on your system.  Note: you can download it [here](http://www.mathworks.com/matlabcentral/fileexchange/31141-felicity).

# Details

Perform the following steps.

* Install a C++ compiler, e.g. MS Visual Studio, or gcc, etc...  See the web site

http://www.mathworks.com/support/compilers/R2012a/win64.html

for more info on installing a compatible C++ compiler with MATLAB.

* Configure the C++ compiler with the command:  `mex -setup` at the MATLAB prompt.
* Unzip the file `FELICITY.zip` to a directory on your computer.
* Open MATLAB and change to the directory you created in the previous step.  Execute the m-script `FELICITY_paths` to add FELICITY to your MATLAB path.
* Execute the script `test_FELICITY` to run a series of unit tests to verify that FELICITY works with your system.
* Execute the script `FELICITY_user_help` to see a listing of FELICITY classes and m-files that are relevant to the user.

# Other Requirements

FELICITY *requires* the MATLAB Symbolic Computing toolbox.

# Troubleshooting

A common problem with using FELICITY is not having the correct C++ compiler linked to MATLAB.  If the compiler is not setup correctly, then FELICITY will not work.  But this is *not* an issue with FELICITY; this is a problem with MATLAB and your compiler.

## Problem On LINUX

For example, this web page describes a common problem people have when running on LINUX:

http://www.mathworks.com/matlabcentral/newsreader/view_thread/162466

If you google this compiler problem, you will find more info.

## Problem On Windows

Compiler problems can also happen with Windows.  For example, several user comments here:

http://www.mathworks.com/matlabcentral/fileexchange/31141-felicity/

describe some issues with installing MS Visual C++ and interfacing with MATLAB.  But again, this is not related to FELICITY.  This is a MATLAB/compiler issue.

## Problem on Mac

Some Mac users have experienced some issues with running FELICITY.  Unfortunately, I don't have a Mac, so I cannot trouble-shoot this.  If anyone can offer advice, please do.

## How To Fix Compiler Problem

There are two ways to fix this:

* Google search for the problem you are having and hopefully find a solution.
* *Complain* to MATLAB support people.  If you are a licensed user, they _must_ help you!  In my experience, they always get back to me in a reasonable time and their advice is useful.

# Other Issues With FELICITY

If you are having problems using the FELICITY toolbox, then either

* send me an email at:   walker@math.lsu.edu; or
* go to the FELICITY forum and post a comment describing the problem:  [TODO]
