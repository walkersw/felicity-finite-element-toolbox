Smoothing Meshes
================

# Introduction

Some types of finite element computations require one to move the mesh, i.e. the positions of vertices of the mesh are changed.  This could be for a time-dependent ALE (Arbitrary-Lagrangian-Eulerian) method.

Because of this, meshes may become distorted with "bad" element qualities.  Usually this means the angles of the elements (triangles or tetrahedrons) are close to being degenerate (near 0 or 180 degrees).  For reasons of accuracy and conditioning, it is best to avoid such degenerate elements.

# Mesh Smoothing by Minimizing E_ODT

FELICITY provides an implementation of a mesh smoother based on the method described in:

Alliez, P.; Cohen-Steiner, D.; Yvinec, M. & Desbrun, M. "Variational Tetrahedral Meshing," ACM Trans. Graph., ACM, 2005, 24, 617-625

and

Chen, L. "Mesh smoothing schemes based on optimal Delaunay triangulations," 13th International Meshing Roundtable, Sandia National Laboratories, 2004, 109-120

# Example: Smooth Deformation

At the MATLAB prompt (or in a script file) type the following:
```matlab
% create mesh of unit square
N_pts = 30+1;
[Elem, Vtx] = bcc_triangle_mesh(N_pts,N_pts);
TR = triangulation(Elem,Vtx);
fB = TR.freeBoundary;
Bdy_Vtx_Indices = unique(fB(:));
Vtx_Attach = TR.vertexAttachments;
clear TR;

disp('Number of Vertices and Elements:');
[size(Vtx,1), size(Elem,1)]
```

Now we apply a smooth deformation to the mesh vertices:
```matlab
% deform mesh
STD = 0.1;
%Gauss = @(x,y) exp(-(x.^2 + y.^2)/STD);
Gauss_x = @(xx) -(2/STD) * (xx(:,1) - 0.5) .* exp(-((xx(:,1) - 0.5).^2 + (xx(:,2) - 0.5).^2)/STD);
Gauss_y = @(xx) -(2/STD) * (xx(:,2) - 0.5) .* exp(-((xx(:,1) - 0.5).^2 + (xx(:,2) - 0.5).^2)/STD);
displace = (0.5*STD) * [Gauss_x(Vtx), Gauss_y(Vtx)];
% set the displacement to zero on the boundary
displace(Bdy_Vtx_Indices,1) = 0;
displace(Bdy_Vtx_Indices,2) = 0;
Vtx_displace = Vtx + displace;
```

To view the mesh, plot it with these commands:
```matlab
figure;
subplot(1,2,1);
trimesh(Elem,Vtx_displace(:,1),Vtx_displace(:,2),0*Vtx_displace(:,2));
view(2);
AX = [0 1 0 1];
axis(AX);
axis equal;
axis(AX);
title('Deformed Mesh');
```

Next, we smooth the mesh by calling the sub-routine `FEL_Mesh_Smooth`:
```matlab
% optimize mesh vertices
Num_Sweeps = 30;
Vtx_Indices_To_Update = setdiff((1:1:size(Vtx_displace,1))',Bdy_Vtx_Indices);
Vtx_Smooth_1 = FEL_Mesh_Smooth(Vtx_displace,Elem,Vtx_Attach,Vtx_Indices_To_Update,Num_Sweeps);
```
Note: calling `FEL_Mesh_Smooth` requires you to give it the mesh connectivity (and initial vertex positions), the elements attached to each vertex (see `Vtx_Attach`), a list of vertex indices to actually move, and the number of times to loop through the vertex list.  In the example above, we did *not* update the vertices on the boundary of the mesh.

Finally, we plot the smoothed (optimized) mesh:
```matlab
subplot(1,2,2);
trimesh(Elem,Vtx_Smooth_1(:,1),Vtx_Smooth_1(:,2),0*Vtx_Smooth_1(:,2));
view(2);
axis(AX);
axis equal;
axis(AX);
title('Smoothed Mesh After 30 Gauss-Seidel Iterations');
```
Notice that it took 30 iterations to achieve a reasonably uniform mesh.  This is because the Gauss-Seidel method is used.

# Example: Random Perturbation

Now lets apply a random displacement of the vertex positions:
```matlab
% perturb mesh
displace = 0.5 * (1/N_pts) * (rand(size(Vtx,1),2) - 0.5);
% set the displacement to zero on the boundary
displace(Bdy_Vtx_Indices,1) = 0;
displace(Bdy_Vtx_Indices,2) = 0;
Vtx_perturb = Vtx + displace;
```
and plot it:
```matlab
figure;
subplot(1,2,1);
trimesh(Elem,Vtx_perturb(:,1),Vtx_perturb(:,2),0*Vtx_perturb(:,2));
view(2);
AX = [0 1 0 1];
axis(AX);
axis equal;
axis(AX);
title('Perturbed Mesh');
```

Next, smooth the mesh:
```matlab
% optimize mesh vertices
Num_Sweeps = 2;
Vtx_Smooth_2 = FEL_Mesh_Smooth(Vtx_perturb,Elem,Vtx_Attach,Vtx_Indices_To_Update,Num_Sweeps);
```
and plot it:
```matlab
subplot(1,2,2);
trimesh(Elem,Vtx_Smooth_2(:,1),Vtx_Smooth_2(:,2),0*Vtx_Smooth_2(:,2));
view(2);
axis(AX);
axis equal;
axis(AX);
title('Smoothed Mesh After 2 Gauss-Seidel Iterations');
```
Notice the smoothed mesh looks very nice even though only two iterations were used.  This is because the Gauss-Seidel method is very effective at removing high frequency components from the solution (which a random perturbation necessarily has).

Note: this method is also implemented for 3-D tetrahedral meshes.  The routine `FEL_Mesh_Smooth` is called in _exactly_ the same way.

This tutorial can be found in ".\Demo\Mesh_Smoothing_2D" sub-directory of FELICITY.