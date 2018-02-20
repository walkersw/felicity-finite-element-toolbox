Mesh Generation and Solving a PDE
=================================

This tutorial takes you from generating a mesh of a non-trivial domain, to solving a PDE (Laplace's equation) on that mesh with non-trivial boundary conditions.

# Laplace's Equation on a Disk With Holes (2D)

Weak Formulation:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Mesh_Gen_With_PDE_Solve_weak_formulation.jpg|width=760|alt=Weak Form Laplace Eqn with Holes]]

Domain \Omega:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Mesh_Gen_With_PDE_Solve_Domain_Fig.jpg|width=670|alt=image of Domain Omega]]

Note: you should do the tutorial [Mesh Generation With TIGER](../wiki/Mesh_Generation_with_TIGER_1) before proceeding.

# Input File

First we need to generate code that will assemble the matrices we need to form the discrete finite element version of the weak formulation described above.

In the MATLAB editor, create the following m-function and name it `MatAssem_Mesh_Gen_With_Solving_PDE.m`: 
```matlab
function MATS = MatAssem_Mesh_Gen_With_Solving_PDE()

% define domain (2-D)
Omega = Domain('triangle');
Gamma = Domain('interval') < Omega;
% there will be other boundaries also, but they are not needed here.

% define finite element spaces
P1_Space = Element(Omega,lagrange_deg1_dim2,1);

% define functions on FE spaces
v = Test(P1_Space);
u = Trial(P1_Space);

% define geometric function on \Gamma
gf = GeoFunc(Gamma);

% define FEM matrices
Stiff_Matrix = Bilinear(P1_Space,P1_Space);
Stiff_Matrix = Stiff_Matrix + Integral(Omega, v.grad' * u.grad );

Neumann_Matrix = Linear(P1_Space);
Neumann_Matrix = Neumann_Matrix + Integral(Gamma, 2*cos(pi*(gf.X(1) + gf.X(2))) * v.val );

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 5;
% define geometry representation - Domain, (default to piecewise linear)
G1 = GeoElement(Omega);
% define a set of matrices
MATS = Matrices(Quadrature_Order,G1);

% collect all of the matrices together
MATS = MATS.Append_Matrix(Stiff_Matrix);
MATS = MATS.Append_Matrix(Neumann_Matrix);

end
```

The Linear form Neumann_Matrix represents the Neumann condition.  Note that we implemented the _g_ function (see PDE formulation above) with
```matlab
2*cos(pi*(gf.X(1) + gf.X(2)))
```
where `gf` is a geometric function on `Gamma` that provides access to geometric information:
* `gf.X(i)` is a symbolic variable that represents the `i`-th coordinate restricted to `Gamma`.  For example, `gf.X(1)` is really the _x_ coordinate evaluated on `Gamma`.
* `gf.T` is a symbolic variable representing the unit tangent vector on `Gamma`.
* `gf.N` is a symbolic variable representing the unit normal vector on `Gamma`.
* etc...
See chapter 4 of the PDF manual for more info.

# Compile It!

Put the file `MatAssem_Mesh_Gen_With_Solving_PDE.m` into a directory that is *in your MATLAB path*.  Compile it by running:
```matlab
Convert_Form_Definition_to_MEX(@MatAssem_Mesh_Gen_With_Solving_PDE,{},'mex_Mesh_Gen_With_Solving_PDE_assemble');
```

Here we named the executable `mex_Mesh_Gen_With_Solving_PDE_assemble`.  See the tutorial [Solve Laplace's Equation](../wiki/Solve_Laplaces_Eqn_1) for more info on this compilation step.

# How To Solve The Problem

We now walk you through creating the mesh, defining the subdomains, defining finite element spaces and boundary conditions, assembling the matrices, solving the system, and plotting the solution.

Note: you can either type the following commands at the MATLAB prompt or put them into a separate script file.

## Create The Mesh

Start by creating a mesh generator:
```matlab
% create the (2-D) mesher object
Box_Dim = [-1, 1];
Num_BCC_Points = 61;
Use_Newton = false;
TOL = 1E-3; % tolerance to use in computing the cut points
% BCC mesh of the square [-1,1] x [-1,1]
MG = Mesher2Dmex(Box_Dim,Num_BCC_Points,Use_Newton,TOL);
```

Now define a level set function that represents our domain \Omega:
```matlab
% create a disk shaped domain with 3 holes in it
LS = LS_Many_Ellipses();
LS.Param.rad_x = [0.7, 0.15, 0.15, 0.15];
LS.Param.rad_y = [0.7, 0.15, 0.15, 0.15];
LS.Param.cx    = [0,  0,    0.35*(sqrt(3)/2), -0.35*(sqrt(3)/2)];
LS.Param.cy    = [0, -0.35, 0.35/2,            0.35/2];
LS.Param.sign  = [1, -1, -1, -1];
```
and create the mesh:
```matlab
% setup up handle to interpolation routine
Interp_Handle = @(pt) LS.Interpolate(pt);

% mesh it!
MG = MG.Get_Cut_Info(Interp_Handle);
[TRI, VTX] = MG.run_mex(Interp_Handle);
```

For convenience, we shall use a FELICITY mesh class to hold the mesh information:
```matlab
% create a Mesh object
Mesh = MeshTriangle(TRI,VTX,'Omega');
% remove any vertices that are not referenced by the triangulation
Mesh = Mesh.Remove_Unused_Vertices;
```
This last line deserves a comment.  The output of the mesh generator often contains vertices that are *not* referenced by the triangulation.  In fact, you may have seen MATLAB give a warning about this.  This is ok.  HOWEVER, you must be careful if you want to use the triangulation data as your DoFmap for a piecewise linear finite element space (which we do below).  In this case, you should use the `Remove_Unused_Vertices` method to delete any unused vertices and renumber the remaining ones; the new mesh object will contain fewer vertices and all of them will be referenced by the triangulation.

Now plot the domain to make sure it looks correct:
```matlab
% plot the domain Omega
figure;
Mesh.Plot;
title('Mesh of Domain \Omega');
AX = 0.8 * [-1 1 -1 1];
axis(AX);
axis equal;
axis(AX);
hold on;
text(0,0.74,'\Gamma');
text(0,-0.25,'\Sigma_1');
text(0.38,0.2,'\Sigma_2');
text(-0.43,0.2,'\Sigma_3');
hold off;
```

## Define Subdomains

Next, we find and store information for the various subdomains.  Note: all of the subdomains here are contained in the boundary of \Omega.

First, we use basic geometric information to identify `Gamma`:
```matlab
% find subdomains
Bdy_of_Omega = Mesh.freeBoundary(); % set of edges
% find the outer part: Gamma
Bdy_of_Omega_MidPt = 0.5 * (Mesh.Points(Bdy_of_Omega(:,1),:) + Mesh.Points(Bdy_of_Omega(:,2),:));
Bdy_of_Omega_MAG = sqrt(sum(Bdy_of_Omega_MidPt.^2,2));
Outer_Mask  = (Bdy_of_Omega_MAG > 0.65);
Gamma_Edges = Bdy_of_Omega(Outer_Mask,:);
```
and we apply similar processing to identify `Sigma_i`:
```matlab
% find the other parts:  Sigma_1, Sigma_2, Sigma_3.
Sigma_1_Mask = (Bdy_of_Omega_MidPt(:,2) < 0) & ~Outer_Mask;
Sigma_2_Mask = (Bdy_of_Omega_MidPt(:,1) > 0) & (Bdy_of_Omega_MidPt(:,2) > 0) & ~Outer_Mask;
Sigma_3_Mask = (Bdy_of_Omega_MidPt(:,1) < 0) & (Bdy_of_Omega_MidPt(:,2) > 0) & ~Outer_Mask;
Sigma_1_Edges = Bdy_of_Omega(Sigma_1_Mask,:);
Sigma_2_Edges = Bdy_of_Omega(Sigma_2_Mask,:);
Sigma_3_Edges = Bdy_of_Omega(Sigma_3_Mask,:);
```

Then, we store this information in the mesh object:
```matlab
% now define subdomains
Mesh = Mesh.Append_Subdomain('1D','Gamma',Gamma_Edges);
Mesh = Mesh.Append_Subdomain('1D','Sigma_1',Sigma_1_Edges);
Mesh = Mesh.Append_Subdomain('1D','Sigma_2',Sigma_2_Edges);
Mesh = Mesh.Append_Subdomain('1D','Sigma_3',Sigma_3_Edges);
```

Lastly, we create embedding data that will be used by the matrix assembler:
```matlab
% create subdomain embedding data
DoI_Names = {'Omega'; 'Gamma'};
Subdomain_Embed = Mesh.Generate_Subdomain_Embedding_Data(DoI_Names);
```

## Define The Finite Element Space

First, define the DoFmap using the triangulation:
```matlab
P1_Space_DoFmap = uint32(Mesh.ConnectivityList);
```
Recall the comment before about `Remove_Unused_Vertices`.

Next, use a FELICITY class to contain information about the finite element space:
```matlab
% define finite element space object
P1_RefElem = ReferenceFiniteElement(lagrange_deg1_dim2());
P1_Lagrange_Space = FiniteElementSpace('Solution', P1_RefElem, Mesh, 'Omega');
% store DoFmap
P1_Lagrange_Space = P1_Lagrange_Space.Set_DoFmap(Mesh,P1_Space_DoFmap);
```
See the tutorial [Managing DoFs: Part 1](../wiki/Managing_DoFs_1) for more info on how `FiniteElementSpace` works.

Then, set subdomains where Dirichlet boundary conditions are set:
```matlab
P1_Lagrange_Space = P1_Lagrange_Space.Set_Fixed_Subdomains(Mesh,{'Sigma_1', 'Sigma_2', 'Sigma_3'});
```
"Fixed" refers to the fact that any Degrees-of-Freedom (DoFs) on those subdomains are fixed (given) by some condition.

## Assemble Matrices

Call the executable that you generated before:
```matlab
FEM = mex_Mesh_Gen_With_Solving_PDE_assemble([],Mesh.Points,P1_Space_DoFmap,[],Subdomain_Embed,P1_Lagrange_Space.DoFmap);
```

## Setup Linear System

Next, we setup the linear system (with boundary conditions) and solve it.

Start with:
```matlab
disp('Setup Laplace equation with Dirichlet conditions on Sigma_1, Sigma_2, Sigma_3:');
Soln = zeros(P1_Lagrange_Space.num_dof,1); % init
A    = FEM(2).MAT;
```
Set the Dirichlet conditions:
```matlab
disp('----> set u = -1 on Sigma_1.');
Sigma_1_DoFs = P1_Lagrange_Space.Get_DoFs_On_Subdomain(Mesh,'Sigma_1');
Soln(Sigma_1_DoFs,1) = -1;
disp('----> set u = +1 on Sigma_2 and Sigma_3.');
Sigma_2_DoFs = P1_Lagrange_Space.Get_DoFs_On_Subdomain(Mesh,'Sigma_2');
Soln(Sigma_2_DoFs,1) = 1;
Sigma_3_DoFs = P1_Lagrange_Space.Get_DoFs_On_Subdomain(Mesh,'Sigma_3');
Soln(Sigma_3_DoFs,1) = 1;
```
Create the right-hand-side data (i.e. the ``load'' vector) to contain the Neumann condition and the Dirichlet conditions:
```matlab
% set right-hand-side (RHS) data
RHS = FEM(1).MAT; % Neumann matrix
RHS = RHS - A * Soln;
```

Get a list of the ``Free'' DoFs:
```matlab
FreeDoFs = P1_Lagrange_Space.Get_Free_DoFs(Mesh);
```
then solve the system just for those DoFs:
```matlab
disp('Solve linear system with backslash:');
Soln(FreeDoFs,1) = A(FreeDoFs,FreeDoFs) \ RHS(FreeDoFs,1);
```

## Plot Solution

Finally, plot the solution:
```matlab
figure;
h1 = trisurf(Mesh.ConnectivityList,Mesh.Points(:,1),Mesh.Points(:,2),Soln);
shading interp;
title('Solution on \Omega','FontSize',14);
set(gca,'FontSize',14);
colorbar;
view(2);
axis(AX);
axis equal;
axis(AX);
```

Note: the solution should look like:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Mesh_Gen_With_PDE_Solve_Solution_Fig.jpg|width=670|alt=image of solution]]