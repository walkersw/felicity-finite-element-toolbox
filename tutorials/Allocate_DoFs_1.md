Allocating Degrees-of-Freedom (DoFs)
====================================

# Introduction

In the previous tutorials, the finite element space was very simple: piecewise linear, continuous Lagrange polynomials.  In this case, the Degrees-of-Freedom (DoFs) correspond directly to the vertices of the mesh.  Hence, you can use the mesh element connectivity data to know which DoFs belong to each element.

This is not true for higher order elements (e.g. piecewise quadratic) or more exotic elements (e.g. Raviart-Thomas).

Allocating DoFs is not just a simple matter of distributing a bunch of indices over a set of mesh elements.  You must ensure that neighboring elements "share" common DoFs.

It is certainly possible for you to write a special purpose MATLAB code to allocate DoFs for a particular element.  However, FELICITY can automatically generate a code to do this for you!

In this tutorial, we will generate code (MEX file) to allocate DoFs for continuous Lagrange, piecewise linear and piecewise quadratic finite element spaces on a 2-D triangle mesh.  Then we will input a sample mesh to this MEX file to test it.

# Code Generation

Type the following at the MATLAB prompt (or put it in a script file):
```matlab
Elem    = lagrange_deg1_dim2(); % piecewise linear
Elem(2) = lagrange_deg2_dim2(); % piecewise quadratic
```

This creates some MATLAB structs that contain important info about the finite element spaces.

Choose a directory on your hard drive to store the MEX file.  Make sure it is in your MATLAB path.  Define a MATLAB string variable that records that directory name:
```matlab
Main_Dir = 'C:\Your_Favorite_Directory\'; 
```

Now generate and compile the DoF allocation code by typing the following command at the MATLAB prompt and press "ENTER":
```matlab
Create_DoF_Allocator(Elem,'mexDoF_Example_2D',Main_Dir);
```

Here we named the executable `mexDoF_Example_2D` (see next section).

# Run It!

First, load up a 2-D triangle mesh:
```matlab
[Vtx, Tri] = Standard_Triangle_Mesh_Test_Data();
tr = triangulation(Tri,Vtx); % triangulation is a built-in MATLAB class
```

For the purposes of illustration, we use the MATLAB class `triangulation` instead of FELICITY's mesh class.

Plot the mesh:
```matlab
figure;
p1 = trimesh(tr.ConnectivityList,tr.Points(:,1),tr.Points(:,2),0*tr.Points(:,2));
view(2);
axis equal;
shading interp;
set(p1,'edgecolor','k'); % make mesh black
title('Sample Mesh for Automatic DoF Allocation Test');
```

Now, allocate the DoFs:
```matlab
[P1_DoFmap, P2_DoFmap] = mexDoF_Example_2D(uint32(tr.ConnectivityList));
```

Note: try running `mexDoF_Example_2D` with the wrong number of arguments; it will give you an error message.  But it will also tell you what the inputs should be, and what the outputs are.

You can now display the DoFmaps:
```matlab
disp('-----------------------');
disp('P1 DoFmap:');
disp(P1_DoFmap);
disp('-----------------------');
disp('P2 DoFmap:');
disp(P2_DoFmap);
disp('-----------------------');
```

`P1_DoFmap` is an Mx3 matrix, where M is the number of triangles.  For the mesh in this example, M = 16.  It has 3 columns because there are only 3 independent basis functions in the piecewise linear polynomial space on each triangle.

`P2_DoFmap` is an Mx6 matrix.  It has 6 columns because there are 6 independent basis functions in the piecewise quadratic polynomial space on each triangle.

Row k (of both matrices) corresponds to triangle k in the mesh.

These DoFmaps can now be used in a matrix assembly routine to generate sparse finite element matrices (see the other tutorials).