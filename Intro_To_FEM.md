Introduction to the Finite Element Method (FEM)
===============================================

# Introduction

The finite element method in its full glory can be very intimidating.  One has to know many things, such as partial differential equation (PDE) theory, numerical analysis, interpolating with basis functions, matrix assembly, solving linear systems, etc...

So here are some suggestions on ways to familiarize yourself with the Finite Element Method (FEM) and how it is implemented.

# Resources Available

## Laplace's Equation in MATLAB

I highly recommend starting by reading this:

J. Alberty, C. Carstensen, S. A. Funken, "Remarks Around 50 Lines Of Matlab: Short Finite Element Implementation," Numerical Algorithms, 1998, 20, 117-137

It can be found here (along with the MATLAB code):

http://www.math.hu-berlin.de/~cc/english/software/shortFE.html

The code solves Laplace's equation in 2-D on unstructured meshes.  The paper also tells you how to modify the code to do more complicated things, like solving non-linear problems and implementing Laplace's equation in 3-D.

This is how I started learning about the finite element method, simultaneous with taking a graduate course on the mathematics of the finite element method.  I also recommend taking a course, or at least some formal self-study.

## Video Lectures on Scientific Computing

For more advanced instruction, with particular emphasis on C++, the following video lectures may be useful to you:

http://www.math.tamu.edu/~bangerth/videos.html

The videos focus on the package deal.ii, but many of the lessons covered are general and pertain to any kind of scientific computing.