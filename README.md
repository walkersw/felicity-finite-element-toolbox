# felicity-finite-element-toolbox

                       The FELICITY Package
 Finite ELement Implementation and Computational Interface Tool for You
                 (C) 01/17/2016, Shawn W. Walker

This code is open source under the BSD license (see the LICENSE.txt file).


DESCRIPTION
========================================================================

This code implements mesh manipulation and refinement, generic finite element matrix assembly for problems in 1D, 2D, and 3D, interpolation of finite element data, as well as some other tools.  See the accompanying FELICITY.pdf and Quickref.pdf files for more info.

USAGE
========================================================================

0. Download the package from:  http://www.mathworks.com/matlabcentral/fileexchange/31141-felicity

1. Extract the given .zip file and read the ``FELICITY.pdf'' file for more info.  In particular, see the INSTALLATION instructions in Chapter 1 of the PDF file.

2. You need a C++ compiler. I use MS Visual C++, express edition.  If you use LINUX, then the gcc compiler should be fine.  You can configure MATLAB to use a compatible C++ compiler by typing ``mex -setup'' at the MATLAB prompt.

3. Check the manual for tutorials.  SWW: still in the process of moving tutorials over to GitHub.

4. Type ``FELICITY_user_help'' at the MATLAB prompt for info on useful classes within FELICITY.
Look in the ``Demo'' directories for examples of how to use FELICITY.  Also, check out the quick refernce guide: Quickref.pdf.

5. Warning!  Make sure that the meshes you use are positively oriented, e.g. make sure edge, triangle, and tetrahedra connectivity lists are positively oriented (i.e. right-hand-rule).  If you do not know what this means, then you should not use this toolbox!


COMPATIBILITY NOTES
========================================================================
The tool was developed in its current form with R2015a.

You need the MATLAB Symbolic Math Toolbox to use FELICITY.
You need a C++ compiler that MATLAB can use with its "mex" command.


Tested on these systems:

-- Windows 7, 10, 64-bit
Fully functional with R2013a, R2015a.

-- LINUX KDE/Ubuntu, 64-bit
Fully functional with R2014a.


BUG REPORTS AND FEEDBACK
========================================================================
Please report any problems and/or bugs to:  walker@math.lsu.edu


ACKNOWLEDGEMENTS
========================================================================

I would like to acknowledge the following contributions from the MATLAB file exchange.
These files are now included in FELICITY.


http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test

by  Darren Engwirda

http://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-volume

by  Sven Holcombe

http://www.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class

by  Oliver Woodford
