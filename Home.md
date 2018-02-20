Welcome to the FELICITY wiki!
=============================

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/FELICITY_Logo_Pic.jpg|alt=Welcome to the FELICITY wiki!]]

FELICITY: Finite ELement Implementation and Computational Interface Tool for You
--------------------------------------------------------------------------------

See the MATLAB central file exchange to download the toolbox (.zip file):

http://www.mathworks.com/matlabcentral/fileexchange/31141-felicity

# Introduction

This is a MATLAB/C++ code for solving PDEs that are discretized by a finite element method on unstructured simplex grids. It uses a Domain-Specific-Language (DSL) to help streamline implementation of FE discretizations (e.g. matrix assembly) by automatic code generation. The resultant sparse matrices can be manipulated by MATLAB for ease in solving a PDE on a triangular (or tetrahedral) mesh.

# Features

* PDEs can be defined on 1-D, 2-D, and 3-D domains. Moreover, the domains can be curves or surfaces embedded in 3-D. For example, you can solve the Laplace-Beltrami equation on a 2-D surface in 3-D.
* Can have multiple interacting sub-domains in a single problem. For example, can have a 1-D curve sub-domain embedded in a 3-D bulk mesh. Can also have a 2-D surface sub-domain embedded in a 3-D bulk-mesh.
* Can define bilinear and linear forms with contributions from multiple embedded sub-domains of different dimensions (i.e. co-dimension >= 0).
* Can do higher order geometry (e.g. quadratic triangle mappings).
* New elements are easily added using a simple ``flat'' m-file.
* Automatically generate custom matrix assembly codes that are callable from MATLAB.
* Automatically generate DoF (Degree-of-Freedom) numbering/allocation for any element implemented in FELICITY.
* Generic mesh classes that implement some useful mesh routines, such as 1-D and 2-D adaptive mesh refinement.
* Some mesh generation utilities (see below).
* MATLAB classes for managing FEM spaces.
* H(div) and H(curl) elements are implemented.
* Support for finite element interpolation in 1-D, 2-D, and 3-D.
* Efficient C++ implementations of bitree, quadtree, and octree; useful for nearest neighbor searching.
* Efficient closest point searching of simplex meshes. This includes finding closest points on surface meshes in 3-D.

Please see the manual (PDF) in the .zip file for more information.

# Citing FELCITY

If you use this toolbox in your work, then you must acknowledge it.  Please cite this wiki website and/or the MATLAB central file exchange URL for FELICITY (bibtex entries are below).  A paper will be forthcoming.
```
@WWW{Walker_FEL_Wiki2017,
  author = {Shawn W. Walker},
  title  = {FELICITY Wiki Documentation},
  year   = {2017},
  url    = {https://github.com/walkersw/felicity-finite-element-toolbox/wiki},
}

@WWW{FELICITY_REF,
  author       = {Shawn W. Walker},
  title        = {{FELICITY}: {F}inite {EL}ement {I}mplementation and {C}omputational {I}nterface {T}ool for {Y}ou},
  url          = {http://www.mathworks.com/matlabcentral/fileexchange/31141-felicity},
  howpublished = {http://www.mathworks.com/matlabcentral/fileexchange/31141-felicity},
}
```

# Hey, I Just Want The Mesh Generator

Then download the FELICITY package and only keep this sub-directory:
```
./FELICITY/Static_Codes/Isosurface_Meshing
```
Next, run `compile_mex_2D_mesh_tiger_code`, `compile_mex_3D_mesh_tiger_code` to compile the C++ code.

Then look at the tutorial: [Mesh Generation With TIGER: Part 1](../wiki/Mesh_Generation_with_TIGER_1).
