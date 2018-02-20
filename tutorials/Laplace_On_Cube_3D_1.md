Solving the 3-D Laplace Equation
================================

This is a tutorial on solving the 3-D Laplace equation with matrix re-assembly, when you have a large number of DoFs.

# Scalar Laplace's Equation on a Cube (3-D)

Weak Formulation:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_On_Cube_3D_weak_formulation.jpg|width=720|alt=Weak Form 3-D Laplace Eqn]]

# Input File

In the MATLAB editor, create the following m-function and name it `MatAssem_Laplace_On_Cube_3D.m`: 

```matlab
function MATS = MatAssem_Laplace_On_Cube_3D()

% define domain (3-D volume)
Omega = Domain('tetrahedron');

% define finite element spaces
P1_Space = Element(Omega,lagrange_deg1_dim3,1);

% define functions on FE spaces
v = Test(P1_Space);
u = Trial(P1_Space);

% define FEM matrices
Stiff_Matrix = Bilinear(P1_Space,P1_Space);
Stiff_Matrix = Stiff_Matrix + Integral(Omega, v.grad' * u.grad );

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 1;
% define geometry representation - Domain, (default to piecewise linear)
G1 = GeoElement(Omega);
% define a set of matrices
MATS = Matrices(Quadrature_Order,G1);

% collect all of the matrices together
MATS = MATS.Append_Matrix(Stiff_Matrix);

end
```

# Compile It!

Put the file *MatAssem_Laplace_On_Cube_3D.m* into a directory that is *in your MATLAB path*.  Compile it by running:
```matlab
Convert_Form_Definition_to_MEX(@MatAssem_Laplace_On_Cube_3D,{},'mex_Laplace_On_Cube_3D_assemble'); 
```

See the tutorial [Solve Laplace's Equation](../wiki/Solve_Laplaces_Eqn_1) for more info on this compilation step.

# Run It!

We will solve the 3-D Laplace problem. First, type the following commands at the MATLAB prompt (or put them into a separate script file):

```matlab
Num_Pts = 20+1;
[Omega_Tet, Omega_Vertex] = regular_tetrahedral_mesh(Num_Pts,Num_Pts,Num_Pts);
Mesh = MeshTetrahedron(Omega_Tet, Omega_Vertex, 'Omega');
clear Omega_Tet Omega_Vertex;
```

Next, define the boundary of \Omega:
```matlab
% define subdomains
Bdy_Faces = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('2D','Bdy',Bdy_Faces);
```

Next, define the local-to-global Degree-of-Freedom map (DoF map) for the finite element space. Since we are using piecewise linear continuous basis functions for the solution, we can simply type:
```matlab
P1_Space_DoFmap = uint32(Mesh.ConnectivityList);
```

Now assemble the matrix from *scratch*:
```matlab
tic
FEM = mex_Laplace_On_Cube_3D_assemble([],Mesh.Points,P1_Space_DoFmap,[],[],P1_Space_DoFmap);
toc
```

Now assemble the matrix *again* using the known sparsity structure:
```matlab
tic
FEM = mex_Laplace_On_Cube_3D_assemble(FEM,Mesh.Points,P1_Space_DoFmap,[],[],P1_Space_DoFmap);
toc
```
Note that we input the`FEM` struct into the first argument of the matrix assembler code.  Also notice that re-assembly is *much faster*.  This is useful when solving *time-dependent* problems.

Next, setup the boundary conditions:
```matlab
Soln = zeros(Mesh.Num_Vtx,1); % init
A    = FEM(1).MAT;
disp('----> set u = sin(2*pi*x) + cos(2*pi*y) + sin(2*pi*z) on Bdy.');
% get Bdy Degrees-of-Freedom (DoF)
Bdy_Nodes = unique(Bdy_Faces(:));
% get Bdy DoF coordinates
Bdy_XC = Mesh.Points(Bdy_Nodes,:);
Soln(Bdy_Nodes,1) = sin(2*pi*Bdy_XC(:,1)) + cos(2*pi*Bdy_XC(:,2)) + sin(2*pi*Bdy_XC(:,3));
RHS = 0*Soln;
RHS = RHS - A * Soln;
```

Finally, solve the system:
```matlab
% get the free nodes of the system
All_Vtx_Indices = (1:1:Mesh.Num_Vtx)';
FreeNodes = setdiff(All_Vtx_Indices,Bdy_Nodes);

disp(' ');
disp(['Size of A = ', num2str(size(A))]);
disp(['Number of Free DoFs = ', num2str(length(FreeNodes))]);
disp(' ');

disp('Solve linear system with backslash:');
tic
Soln(FreeNodes,1) = A(FreeNodes,FreeNodes) \ RHS(FreeNodes,1);
toc
```

NOTE: solving 3-D problems with backslash is rather limiting.  Depending on your machine, you will not be able to make the system matrix too large without running out of memory when executing "backslash".  Moreover, even if it does work, backslash is a *slow* solver for large 3-D problems.  You need an iterative solver.

One option is to use AGMG:

http://homepages.ulb.ac.be/~ynotay/AGMG/

It has a MATLAB interface.  If you have AGMG successfully installed on your system (with the MATLAB interface) then you can do this:
```matlab
% Note: you must have AGMG installed to run this!
disp(' ');
TOL = 1e-8;
disp(['Solve linear system with AGMG with TOL = ', num2str(TOL,'%1.2G'), ':']);
Soln_AGMG = Soln; % keep boundary conditions
tic
Soln_AGMG(FreeNodes,1) = agmg(A(FreeNodes,FreeNodes),RHS(FreeNodes,1),[],TOL);
toc

Soln_Diff = max(abs(Soln - Soln_AGMG));
disp(' ');
disp('Max difference between backslash and AGMG:');
Soln_Diff
```
For a large 3-D elliptic PDE problem, AGMG is *much faster* than backslash.  Note: AGMG is *not* included with FELICITY; you must download it and install it separately.

# Plot It!

You can plot the solution with the following commands:
```matlab
figure;
h1 = trisurf(Bdy_Faces,Mesh.Points(:,1),Mesh.Points(:,2),Mesh.Points(:,3),Soln(:,1));
title('Boundary Mesh and Boundary Condition','FontSize',12);
set(gca,'FontSize',14);
AX = [0 1 0 1 0 1];
axis(AX);
axis equal;
axis(AX);
```

This tutorial can be found in ".\Demo\Laplace_On_Cube_3D" sub-directory of FELICITY.