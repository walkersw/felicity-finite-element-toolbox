Solve the Stokes Equations in 2-D
=================================

# Weak Formulation:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Demo_Stokes_2D_weak_formulation.jpg|width=800|alt=Weak Form 2-D Stokes Eqn]]

We will use a mixed finite element method to approximate the solution to this problem, i.e. we have two finite element spaces of different degree.  One for the velocity field and one for the pressure field.

The pressure space will be approximated by a scalar piecewise linear continuous Lagrange space; we did something similar in [Solve Simple Elasticity Model](../wiki/Solve_Simple_Elasticity_3D_1) for the 3-D case.

However, because pressure is approximated by degree 1 polynomials, we must approximate velocity by degree 2 polynomials (piecewise quadratic).  This is to satisfy the LBB (stability) condition.

# Allocate DoFs for Piecewise Quadratic Polynomials

You should do the tutorial [Allocating Degrees-of-Freedom (DoFs)](../wiki/Allocate_DoFs_1) before proceeding.

Allocating the Degrees-of-Freedom (DoFs) for the piecewise quadratic space is not as easy as for the piecewise linear case.  This is because some DoFs are associated with vertices, and others are associated with edges of the mesh.  You can certainly write a MATLAB code to do this, but it is simpler to just let FELICITY do it for you.

Type the following at the MATLAB prompt (or put it in a script file):
```matlab
Main_Dir = 'C:\Your_Favorite_Directory\'; 
Elem = lagrange_deg2_dim2(); % piecewise quadratic
Create_DoF_Allocator(Elem,'mex_DoF_Lagrange_P2_Allocator_2D',Main_Dir);
```

Here we named the executable `mex_DoF_Lagrange_P2_Allocator_2D`.  See [Allocating Degrees-of-Freedom (DoFs)](../wiki/Allocate_DoFs_1) for more details.

# Input File

In the MATLAB editor, create the following m-function and name it `MatAssem_Stokes_2D.m`: 

```matlab
function MATS = MatAssem_Stokes_2D()

% define domain (2-D)
Omega = Domain('triangle'); % fluid domain (unit square)
Outlet = Domain('interval') < Omega; % right side of \Omega

% define finite element spaces
Vector_P2 = Element(Omega,lagrange_deg2_dim2,2); % 2-D vector valued
Scalar_P1 = Element(Omega,lagrange_deg1_dim2,1);

% define functions on FE spaces
vv = Test(Vector_P2);
uu = Trial(Vector_P2);

q = Test(Scalar_P1);
p = Trial(Scalar_P1);

BC_Out = Coef(Vector_P2);

% define FEM matrices
BC_Matrix = Linear(Vector_P2);
BC_Matrix = BC_Matrix + Integral(Outlet, BC_Out.val' * vv.val );

Div_Pressure = Bilinear(Scalar_P1,Vector_P2);
Div_Pressure = Div_Pressure + Integral(Omega,-q.val * (uu.grad(1,1) + uu.grad(2,2)));

Stress_Matrix = Bilinear(Vector_P2,Vector_P2);
% symmetrized gradient
D_u = uu.grad + uu.grad';
Stress_Matrix = Stress_Matrix + Integral(Omega,sum(sum(D_u .* vv.grad)));

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 4;
% define geometry representation - Domain, (default to piecewise linear)
G1 = GeoElement(Omega);
% define a set of matrices
MATS = Matrices(Quadrature_Order,G1);

% collect all of the matrices together
MATS = MATS.Append_Matrix(BC_Matrix);
MATS = MATS.Append_Matrix(Div_Pressure);
MATS = MATS.Append_Matrix(Stress_Matrix);

end
```

The space `Vector_P2` is for the velocity field.  `Scalar_P1` is for the pressure.

The Linear form `BC_Matrix` represents the stress boundary condition.  The `Stress_Matrix` is a symmetric Bilinear form and is the first diagonal block of the saddle point system.  `Div_Pressure` is a non-symmetric Bilinear form and is the off-diagonal block that captures the divergence free condition.

`BC_Out` is a coefficient function that will contain the stress boundary condition.

# Compile It!

Put the file `MatAssem_Stokes_2D.m` into a directory that is *in your MATLAB path*.  Compile it by running:
```matlab
Convert_Form_Definition_to_MEX(@MatAssem_Stokes_2D,{},'mex_Stokes_2D_assemble'); 
```

Here we named the executable `mex_Stokes_2D_assemble`.  See the tutorial [Solve Laplace's Equation](../wiki/Solve_Laplaces_Eqn_1) for more info on this compilation step.

# Run It!

We will solve the 2-D Stokes problem. First, type the following commands at the MATLAB prompt (or put them into a separate script file):

```matlab
% define square mesh [0, 1] x [0, 1]
[Omega_Tri, Omega_Vertex] = bcc_triangle_mesh(5,5);
Mesh = MeshTriangle(Omega_Tri, Omega_Vertex, 'Omega');
clear Omega_Tri Omega_Vertex;
Mesh = Mesh.Refine;
% plot it
Mesh.Plot;
```

Next, define various sub-domains of \Omega:
```matlab
% define \partial \Omega
Bdy_Edges = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('1D','Bdy',Bdy_Edges);

% get centroids of boundary edges
EE_center = (1/2) * (Mesh.Points(Bdy_Edges(:,1),:) + Mesh.Points(Bdy_Edges(:,2),:));

% find edges on top
Top_Mask = (EE_center(:,2) > 1 - 1e-5);
Top_Edges = Bdy_Edges(Top_Mask,:);
Mesh = Mesh.Append_Subdomain('1D','Top',Top_Edges);
% find edges on bottom
Bot_Mask = (EE_center(:,2) < 0 + 1e-5);
Bot_Edges = Bdy_Edges(Bot_Mask,:);
Mesh = Mesh.Append_Subdomain('1D','Bottom',Bot_Edges);
% find edges on inlet
In_Mask = (EE_center(:,1) < 0 + 1e-5);
In_Edges = Bdy_Edges(In_Mask,:);
Mesh = Mesh.Append_Subdomain('1D','Inlet',In_Edges);
% find edges on outlet
Out_Mask = (EE_center(:,1) > 1 - 1e-5);
Out_Edges = Bdy_Edges(Out_Mask,:);
Mesh = Mesh.Append_Subdomain('1D','Outlet',Out_Edges);
```

Next, we need the local-to-global Degree-of-Freedom map (DoF map) for the finite element spaces. Since we are using piecewise linear continuous basis functions for the pressure, we can simply type:
```matlab
Pressure_DoFmap = uint32(Mesh.ConnectivityList);
```

For the velocity space, it will behoove us to use another FELICITY class `FiniteElementSpace`.  Enter the following into MATLAB:
```matlab
P2_RefElem = ReferenceFiniteElement(lagrange_deg2_dim2());
P2_Lagrange_Space = FiniteElementSpace('Velocity', P2_RefElem, Mesh, 'Omega', 2); % 2 components
```
This class allows us to extract DoFs that are attached to embedded sub-domains.

Now, allocate the DoFs for the velocity field and store them in the `P2_Lagrange_Space` object:
```matlab
Lag_P2_DoFmap = uint32(mex_DoF_Lagrange_P2_Allocator_2D(uint32(Mesh.ConnectivityList)));
P2_Lagrange_Space = P2_Lagrange_Space.Set_DoFmap(Mesh,Lag_P2_DoFmap);
```
Note: only the first component (x-component) of the velocity Degrees-of-Freedom (DoFs) are explicitly stored.  Obtaining the y-component DoFs is easy: just shift the indices by the number of x-component nodes (see the PDE solve part of this tutorial).

Create sub-domain embedding data structures.  This is so the matrix assembly code knows how Outlet is embedded in \Omega:
```matlab
DoI_Names = {'Omega'; 'Outlet'}; % domains of integration
Subdomain_Embed = Mesh.Generate_Subdomain_Embedding_Data(DoI_Names);
```

Next, define the stress boundary condition:
```matlab
% get P2_Lagrange Node coordinates
P2_X = P2_Lagrange_Space.Get_DoF_Coord(Mesh);
% set stress boundary condition
BC_Out_Values = 0*P2_X; % set to the zero vector
% use the following if you want the stress boundary condition
%     to be  (0, sin(2*pi*y))
%BC_Out_Values(:,2) = sin(2*pi*P2_X(:,2));
```

Now assemble the matrices:
```matlab
FEM = mex_Stokes_2D_assemble([],Mesh.Points,uint32(Mesh.ConnectivityList),[],Subdomain_Embed,...
      Pressure_DoFmap,P2_Lagrange_Space.DoFmap,BC_Out_Values);
```
and put the `FEM` data into a special object to allow for easy access of the sparse matrix data:
```matlab
Stokes_Matrices = FEMatrixAccessor('Stokes',FEM);
```

Finally, we obtain the solution of the Stokes system by the following commands.

Get DoFs:
```matlab
Vel_DoFs      = unique(P2_Lagrange_Space.DoFmap(:));
Pressure_DoFs = unique(Pressure_DoFmap(:));
Num_Scalar_Vel_DoF = length(Vel_DoFs);
Num_Pressure_DoF   = length(Pressure_DoFs);
```
Note: `Vel_DoFs` refers to the x-component of the velocity.  The y-component DoFs are obtained by shifting:
```matlab
Vel_DoFs + Num_Scalar_Vel_DoF
```

Initialize solution arrays:
```matlab
Vel_Soln      = zeros(Num_Scalar_Vel_DoF,2); % 2 columns <-> 2-D vector
Pressure_Soln = zeros(Num_Pressure_DoF,1);
```
Access finite element matrices:
```matlab
A = Stokes_Matrices.Get_Matrix('Stress_Matrix');
B = Stokes_Matrices.Get_Matrix('Div_Pressure');
```
Note: by using `Stokes_Matrices`, we do not need to know the *order* of the matrices in the `FEM` struct.  This is especially useful if you must modify the input file `MatAssem_Stokes_2D` to define *additional* bilinear or linear forms (matrices).

Now setup the system matrix:
```matlab
% create saddle point system
Saddle = [A, B'; B, sparse(Num_Pressure_DoF,Num_Pressure_DoF)];
```
Note: because of the way that we ordered the matrices here, the global pressure DoFs are actually
```matlab
Pressure_DoFs + 2*Num_Scalar_Vel_DoF
```

Set the boundary conditions:
```matlab
% compute the parabolic profile over all velocity DoFs
Parabolic_Profile = P2_X(:,2).*(1 - P2_X(:,2)); % y * (1 - y)
% get the DoFs associated with the Inlet
Inlet_DoFs = P2_Lagrange_Space.Get_DoFs_On_Subdomain(Mesh,'Inlet');
% set the x-component of the velocity solution
Vel_Soln(Inlet_DoFs,1) = Parabolic_Profile(Inlet_DoFs,1);
% zero velocity on top and bottom is already set

% put the velocity and pressure solution into a single column vector
Soln = [Vel_Soln(:); Pressure_Soln];
```

Set the right-hand-side data for the matrix system
```matlab
% concatenate stress boundary condition and zero divergence data
R1 = Stokes_Matrices.Get_Matrix('BC_Matrix');
RHS = [R1; zeros(Num_Pressure_DoF,1)];
% include boundary conditions
RHS = RHS - Saddle * Soln;
```

Eliminate the fixed DoFs of the system:
```matlab
Top_DoFs = P2_Lagrange_Space.Get_DoFs_On_Subdomain(Mesh,'Top');
Bot_DoFs = P2_Lagrange_Space.Get_DoFs_On_Subdomain(Mesh,'Bottom');
% use unique b/c some of the nodes may be shared between Top and Inlet for example
Fixed_Scalar_Nodes = unique([Inlet_DoFs; Top_DoFs; Bot_DoFs]);
% get the global velocity nodes that are fixed (Dirichlet condition)
Fixed_Vector_Nodes = [Fixed_Scalar_Nodes; Fixed_Scalar_Nodes + Num_Scalar_Vel_DoF];
% get the remaining free nodes of the system
FreeNodes = setdiff((1:1:length(Soln))', Fixed_Vector_Nodes);
```

Now solve it:
```matlab
% Use backslash to solve
Soln(FreeNodes,1) = Saddle(FreeNodes,FreeNodes) \ RHS(FreeNodes,1);
% now parse the solution back into the velocity and pressure variables
Vel_Soln(:)   = Soln(1:2*Num_Scalar_Vel_DoF,1);
Pressure_Soln = Soln(2*Num_Scalar_Vel_DoF+1:end,1);
```

# Plot It

```matlab
figure;
subplot(2,2,1);
h1 = Mesh.Plot;
title('Omega','FontSize',14);
set(gca,'FontSize',14);
AX = [0 1 0 1];
axis(AX);
axis equal;
axis(AX);

subplot(2,2,2);
h2 = quiver(P2_X(:,1),P2_X(:,2),Vel_Soln(:,1),Vel_Soln(:,2));
title('Vector Velocity Solution','FontSize',14);
set(gca,'FontSize',14);
axis(AX);
axis equal;
axis(AX);

subplot(2,2,3);
h3 = trisurf(Mesh.ConnectivityList,Mesh.Points(:,1),Mesh.Points(:,2),Pressure_Soln);
shading interp;
title('Pressure Solution','FontSize',14);
set(gca,'FontSize',14);
colorbar;
axis(AX);
axis equal;
axis(AX);
```

This tutorial can be found in ".\Demo\Stokes_2D" sub-directory of FELICITY.