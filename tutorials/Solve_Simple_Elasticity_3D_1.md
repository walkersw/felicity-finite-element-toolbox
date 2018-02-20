Solve Laplace's Equation in 3-D
===============================

# Weak Formulation:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Simple_Elasticity_3D_weak_formulation.jpg|width=700|alt=Weak Form 3-D Laplace]]

Note: this is a simple model of elasticity in 3-D.

# Input File

In the MATLAB editor, create the following m-function and name it `MatAssem_Simple_Elasticity_3D.m`:

```matlab
function MATS = MatAssem_Simple_Elasticity_3D()

% define domain (3-D volume)
Omega = Domain('tetrahedron'); % ``elastic'' domain
Sigma = Domain('interval') < Omega; % edge of \Omega

% define finite element spaces
Vector_P1 = Element(Omega,lagrange_deg1_dim3,3); % 3-D vector valued

% define functions on FE spaces
vv = Test(Vector_P1);
uu = Trial(Vector_P1);

Displace = Coef(Vector_P1);

% define geometric function on \Sigma
gSigma = GeoFunc(Sigma);

% define FEM matrices
Stiff_Matrix = Bilinear(Vector_P1,Vector_P1);
Stiff_Matrix = Stiff_Matrix + Integral(Omega,sum(sum(vv.grad .* uu.grad)));

Eval_Matrix = Real(1,1);
Eval_Matrix = Eval_Matrix + Integral(Sigma, Displace.val' * gSigma.T );

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 3;
% define geometry representation - Domain, (default to piecewise linear)
G1 = GeoElement(Omega);

% collect all of the matrices together
MATS = Matrices(Quadrature_Order,G1);
MATS = MATS.Append_Matrix(Eval_Matrix);
MATS = MATS.Append_Matrix(Stiff_Matrix);

end 
```

Here we define the main bilinear form and we define the integral quantity in the "matrix" `Eval_Matrix`.  Note that the finite element space is vector-valued.  We also introduced a geometric function object (for the sub-domain \Sigma) gSigma.  The tangent vector of \Sigma is accessed by gSigma.T.

# Compile It!

Put the file `MatAssem_Simple_Elasticity_3D.m` into a directory that is *in your MATLAB path*.  Now compile it by typing the following command at the MATLAB prompt and press "ENTER":

```matlab
Convert_Form_Definition_to_MEX(@MatAssem_Simple_Elasticity_3D,{},'Assemble_Simple_Elasticity_3D'); 
```

Here we named the executable `Assemble_Simple_Elasticity_3D` (see next section).

# Run It!

We will solve the 3-D Laplace problem. First, type the following commands at the MATLAB prompt (or put them into a separate script file):

```matlab
% create mesh for elastic column
[Omega_Tet, Omega_Vertex] = regular_tetrahedral_mesh(2+1,2+1,10+1);
% scale it
Omega_Vertex(:,1:2) = (1/5) * Omega_Vertex(:,1:2);
Mesh = MeshTetrahedron(Omega_Tet, Omega_Vertex, 'Omega');
% plot it
Mesh.Plot;
```

This defines a tetrahedral mesh for \Omega.

Next, define various sub-domains of \Omega:
```matlab
% \partial \Omega
Bdy_Faces = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('2D','Bdy',Bdy_Faces);
% get centroids of boundary triangles
FF_center = (1/3) * (Mesh.Points(Bdy_Faces(:,1),:) + Mesh.Points(Bdy_Faces(:,2),:) + Mesh.Points(Bdy_Faces(:,3),:));
% find faces on top
Top_Mask = (FF_center(:,3) > 1 - 1e-5);
Top_Faces = Bdy_Faces(Top_Mask,:);
Mesh = Mesh.Append_Subdomain('2D','Top',Top_Faces);
% find faces on bottom
Bot_Mask = (FF_center(:,3) < 0 + 1e-5);
Bot_Faces = Bdy_Faces(Bot_Mask,:);
Mesh = Mesh.Append_Subdomain('2D','Bottom',Bot_Faces);
```

Also, define the sub-domain \Sigma:
```matlab
% find vertices on the edge segment: (0,0,1) - (0,0,0)
Col_Edge_Vtx_Mask = (Mesh.Points(:,1) < 0 + 1e-5) & (Mesh.Points(:,2) < 0 + 1e-5);
Col_Edge_Vtx = Mesh.Points(Col_Edge_Vtx_Mask,:);
All_Vtx_Indices = (1:1:Mesh.Num_Vtx)';
Col_Edge_Vtx_Indices = All_Vtx_Indices(Col_Edge_Vtx_Mask);
[Col_Edge_Vtx, I1] = sortrows(Col_Edge_Vtx,3); % sort in ascending order along z
Col_Edge_Vtx_Indices = Col_Edge_Vtx_Indices(I1);
% define the column edges
Col_Edge_Data = [Col_Edge_Vtx_Indices(1:end-1,1), Col_Edge_Vtx_Indices(2:end,1)];
Mesh = Mesh.Append_Subdomain('1D','Sigma',Col_Edge_Data);
```

Next, we need the local-to-global Degree-of-Freedom map (DoF map) for the finite element space. Since we are using piecewise linear continuous basis functions, we can simply type:
```matlab
P1_DoFmap = uint32(Mesh.ConnectivityList);
```

Define the values of a vector-valued `Displace' coefficient function:
```matlab
Displace = 0*Mesh.Points; % the zero function
```

Create sub-domain embedding data structures.  This is so the matrix assembly code knows how \Sigma is embedded in \Omega:
```matlab
DoI_Names = {'Omega'; 'Sigma'}; % domains of integration
Subdomain_Embed = Mesh.Generate_Subdomain_Embedding_Data(DoI_Names);
```

Now assemble the matrices:

```matlab
FEM = Assemble_Simple_Elasticity_3D([],Mesh.Points,uint32(Mesh.ConnectivityList),[],Subdomain_Embed,P1_DoFmap,Displace);
```

`FEM` is a MATLAB struct that contains the sparse matrices. Finally, we obtain the solution of the PDE by the following commands.

Initialize the solution and set boundary conditions:
```matlab
Soln = zeros(Mesh.Num_Vtx,3); % init the solution
A   = FEM(2).MAT;
% Degree-of-Freedom associated with the Top sub-domain
Top_Nodes = unique(Top_Faces(:));
% set Dirichlet conditions
Soln(Top_Nodes,1) = 0.2;
Soln(Top_Nodes,2) = 0.2;
Soln(Top_Nodes,3) = 0.1;
Soln = Soln(:); % collapse to a single column vector
RHS = 0*Soln;
RHS = RHS - A * Soln; % modify right hand side
% get bottom DoF
Bot_Nodes = unique(Bot_Faces(:));
% bottom nodes are already set to zero

% get the free nodes of the system (for each component)
FreeNodes_X = setdiff(All_Vtx_Indices,[Top_Nodes; Bot_Nodes]);
FreeNodes_Y = [FreeNodes_X + Mesh.Num_Vtx];
FreeNodes_Z = [FreeNodes_Y + Mesh.Num_Vtx];
% combine them all into one array
FreeNodes = [FreeNodes_X; FreeNodes_Y; FreeNodes_Z];
```

Now solve it:
```matlab
Soln(FreeNodes,1) = A(FreeNodes,FreeNodes) \ RHS(FreeNodes,1);
Displace = zeros(Mesh.Num_Vtx,3);
Displace(:) = Soln; % put displacement solution back into 3-D vector form
```

Next, compute the integral quantity
```matlab
Recalc_FEM = Assemble_Simple_Elasticity_3D([],Mesh.Points,uint32(Mesh.ConnectivityList),...
             [],Subdomain_Embed,P1_DoFmap,Displace);

% value of this 1x1 matrix should be 0.05
Integral_Quantity = Recalc_FEM(1).MAT;
```

# Plot It

```matlab
figure;
subplot(1,2,1);
h1 = Mesh.Plot;
title('Reference Mesh','FontSize',14);
set(gca,'FontSize',14);
AX = [0 0.4 0 0.4 0 1.2];
axis(AX);
axis equal;
axis(AX);

subplot(1,2,2);
New_Mesh = Mesh.Set_Points(Mesh.Points + Displace); % make a new displaced mesh
h2 = New_Mesh.Plot;
title('Deformed Mesh','FontSize',14);
set(gca,'FontSize',14);
axis(AX);
axis equal;
axis(AX);
```

For more information, see the PDF manual and the ".\Demo\Simple_Elasticity_3D\" sub-directory of FELICITY.