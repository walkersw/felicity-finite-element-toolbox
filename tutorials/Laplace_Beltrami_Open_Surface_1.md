Solving Laplace-Beltrami on a Higher Order Mesh
===============================================

This tutorial describes how to solve the Laplace-Beltrami equation on a 3-D surface with open boundary using higher order meshes (e.g. piecewise quadratic triangles) instead of the standard piecewise linear case.  This problem is described in the FELICITY paper (coming soon).  Also see the chapter on "Helper Routines For Finite Element Codes" in the FELICITY manual for details on higher order meshes.

Note: this example can be run using the script `test_Laplace_Beltrami_Open_Surface.m` located in
```
./FELICITY/Demo/Laplace_Beltrami_Open_Surface/
```

# Problem Definition

The strong formulation of the PDE is:
[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Strong_Form.jpg|width=320|alt=Strong Formulation]]

Here is an illustration of the domain:
<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Diagram.jpg" width="550" alt="Saddle Domain Diagram">

The *weak* formulation of the PDE is:
[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Weak_Form.jpg|width=320|alt=Weak Formulation]]

# Domain Definition

We will define a piecewise quadratic triangulation of the following parameterized domain:
[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Param.jpg|width=320|alt=Saddle Domain Param]]

See the figure above.

# Curved Finite Element Spaces

We will define a continuous piecewise quadratic (Lagrange) finite element space on the piecewise quadratic mesh, a so-called iso-parametric FE space.  In addition, we will define a piecewise linear space on the piecewise quadratic mesh.  The exact mathematical formulation is:
[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_FE_Spaces.jpg|width=320|alt=Curved FE Spaces]]

# Finite Element Formulation

We wish to solve the following problem:
<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Discrete_Form.jpg" width="600" alt="Discrete Formulation">

## Finite Element Forms

We need the following bilinear forms to construct the linear system representing the finite element formulation above.

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Bilinear.jpg" width="350" alt="Bilinear Forms">

## Abstract Definition File

The (abstract) m-file that defines the bilinear forms is:
```
function MATS = MatAssem_Laplace_Beltrami_Open_Surface()

% define domain (2-D surface in 3-D)
Gamma = Domain('triangle',3);
dGamma = Domain('interval') < Gamma; % subset

% define finite element spaces
V_h = Element(Gamma, lagrange_deg2_dim2,1); % piecewise quadratic
W_h = Element(dGamma,lagrange_deg1_dim1,1); % piecewise linear

% define functions in FE spaces
v_h = Test(V_h);
u_h = Trial(V_h);
mu_h = Test(W_h);

% define (discrete) forms
M = Bilinear(V_h,V_h);
M = M + Integral(Gamma, v_h.val * u_h.val );

K = Bilinear(V_h,V_h);
K = K + Integral(Gamma, v_h.grad' * u_h.grad );

B = Bilinear(W_h,V_h);
B = B + Integral(dGamma, mu_h.val * u_h.val );

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 10;
% define geometry representation - Domain, piecewise quadratic representation
G1 = GeoElement(Gamma,lagrange_deg2_dim2);
% define a set of matrices
MATS = Matrices(Quadrature_Order,G1);

% collect all of the matrices together
MATS = MATS.Append_Matrix(B);
MATS = MATS.Append_Matrix(K);
MATS = MATS.Append_Matrix(M);

end
```
Note that the space `V_h` is piecewise quadratic, and the `GeoElement` is defined to be piecewise quadratic. `W_h` is piecewise linear.

# Convergence Check

We will check the accuracy of our finite element method via a convergence check.

## Numerical Errors

We will evaluate the following quantities of interest in order to compute the convergence rate of the method:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Errors.jpg" width="700" alt="Numerical Errors">

## Abstract Definition File

The (abstract) m-file that defines the numerical errors is:
```
function MATS = Compute_Errors_Laplace_Beltrami_Open_Surface()

% define domain (2-D surface in 3-D)
Gamma = Domain('triangle',3);
dGamma = Domain('interval') < Gamma; % subset

% define finite element spaces
V_h = Element(Gamma, lagrange_deg2_dim2,1); % piecewise quadratic
W_h = Element(dGamma,lagrange_deg1_dim1,1); % piecewise linear
G_h = Element(Gamma, lagrange_deg2_dim2,3); % vector piecewise quadratic

% define functions in FE spaces
u_h = Coef(V_h);
lambda_h = Coef(W_h);

NV_h = Coef(G_h);
TV_h = Coef(G_h);

% geometry access functions
gf_Gamma = GeoFunc(Gamma);
gf_dGamma = GeoFunc(dGamma);

% define exact solns
u_exact = @(u,v) cos(v.*pi.*2.0).*sin(u.*pi.*2.0);
lambda_exact = <too long to show>;

% define (discrete) forms
u_L2_Error_sq = Real(1,1);
u_L2_Error_sq = u_L2_Error_sq + Integral(Gamma, (u_h.val - u_exact(gf_Gamma.X(1),gf_Gamma.X(2)))^2 );

lambda_L2_Error_sq = Real(1,1);
lambda_L2_Error_sq = lambda_L2_Error_sq + Integral(dGamma, (lambda_h.val - lambda_exact(gf_dGamma.X(1),gf_dGamma.X(2)))^2 );

% check normal vector
NV_Error_sq = Real(1,1);
NV_Error_sq = NV_Error_sq + Integral(dGamma, sum((gf_Gamma.N - NV_h.val).^2) );

% check tangent vector
TV_Error_sq = Real(1,1);
TV_Error_sq = TV_Error_sq + Integral(dGamma, sum((gf_dGamma.T - TV_h.val).^2) );

% check Neumann data
xi_h = cross(gf_dGamma.T,gf_Gamma.N);
Neumann_L2_Error_sq = Real(1,1);
Neumann_L2_Error_sq = Neumann_L2_Error_sq + Integral(dGamma, ( lambda_h.val - ( -dot(xi_h,u_h.grad) ) )^2 );

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 10;
% define geometry representation - Domain, piecewise quadratic representation
G1 = GeoElement(Gamma,lagrange_deg2_dim2);
% define a set of matrices
MATS = Matrices(Quadrature_Order,G1);

% collect all of the matrices together
MATS = MATS.Append_Matrix(u_L2_Error_sq);
MATS = MATS.Append_Matrix(lambda_L2_Error_sq);
MATS = MATS.Append_Matrix(NV_Error_sq);
MATS = MATS.Append_Matrix(TV_Error_sq);
MATS = MATS.Append_Matrix(Neumann_L2_Error_sq);

end
```
Note that the space `G_h` corresponds to the representation of the piecewise quadratic triangulation of the domain.

# Implementation

After compiling the abstract definition files (for example, see [Solve Laplace's Equation](../wiki/Solve_Laplaces_Eqn_1)), we now implement a script for assembling the previously defined forms (FE matrices) and quantities of interest.  This involves doing a convergence check, which means we will loop over several meshes of decreasing mesh size.

Note: this can be found in `Execute_Laplace_Beltrami_Open_Surface.m`.

## Initialization

We start with some initialization:
```
Main_Dir = <choose a dir in your MATLAB path>;

Refine_Vec = [2 3 4 5 6 7];
Num_Refine = length(Refine_Vec);

% init arrays to store numerical errors
u_Error_Vec = zeros(Num_Refine,1);
lambda_Error_Vec = zeros(Num_Refine,1);
NV_Error_Vec = zeros(Num_Refine,1);
TV_Error_Vec = zeros(Num_Refine,1);
Neumann_Error_Vec = zeros(Num_Refine,1);
```

Now define the domain parameterization:
```
% define domain parameterization
Psi = @(q1,q2) [q1, q2, 0.5*(q1.^2 - q2.^2)];
```

Next, start the for loop:
```
for kk = 1:Num_Refine
```

## Create The Base (Piecewise Linear) Mesh

Now create a mesh of the unit disk and map it with the parameterization:
```
% BEGIN: define saddle surface (piecewise linear mesh)
Refine_Level = Refine_Vec(kk);
% create mesh of unit disk centered at origin
[Tri, Pts_q] = triangle_mesh_of_disk([0 0 0],1,Refine_Level);

% apply parameterization
Pts_x_y_z = Psi(Pts_q(:,1), Pts_q(:,2));
Mesh = MeshTriangle(Tri, Pts_x_y_z, 'Gamma');
% END: define saddle surface (piecewise linear mesh)
```
This gives us a base piecewise *linear* mesh.

Add the boundary `dGamma` as a sub-domain:
```
% add the open boundary
BDY  = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('1D','dGamma',BDY);
```

## Define a GeoElementSpace

Next, we create a special object for handling the higher order mesh geometry:
```
% create GeoElementSpace
P2_RefElem = ReferenceFiniteElement(lagrange_deg2_dim2());
G_Space = GeoElementSpace('G_h',P2_RefElem,Mesh);
G_DoFmap = DEMO_LB_mex_V_Space_DoF_Allocator(uint32(Mesh.ConnectivityList));
G_Space = G_Space.Set_DoFmap(Mesh,uint32(G_DoFmap));
```
Note: at this point, we have already created a mex file to allocate DoFs for the piecewise quadratic space.

## Create the Piecewise Quadratic Mapping

With this, we get an initial piecewise quadratic FE function that represents the identity map for our mesh:
```
% BEGIN: define the higher order surface
Geo_Points_hat = G_Space.Get_Mapping_For_Piecewise_Linear_Mesh(Mesh);
% initially starts as a piecewise quadratic polynomial (P2 function) interpolating a
% piecewise linear (P1) function (representing the PW linear saddle surface mesh)
```
Since the base mesh is piecewise *linear*, this "quadratic" function is actually interpolating a piecewise linear function.

Now, we create a new variable `Geo_Points` that represents the true piecewise quadratic mapping that we need to define our curved domain:
```
% map back to (q_1,q_2) plane
Geo_q = Geo_Points_hat(:,1:2);
% Geo_q is a piecewise linear approximation of the unit disk
dGamma_DoF_Indices = G_Space.Get_DoFs_On_Subdomain(Mesh,'dGamma');
Geo_q_on_bdy = Geo_q(dGamma_DoF_Indices,:);
% ensure the (q_1,q_2) points on the bdy of the mesh lie *exactly* on the
% boundary of the unit disk: {q_1^2 + q_2^2 < 1}
Norm_on_bdy = sqrt(sum(Geo_q_on_bdy.^2,2));
Geo_q_on_bdy(:,1) = Geo_q_on_bdy(:,1) ./ Norm_on_bdy; % normalize
Geo_q_on_bdy(:,2) = Geo_q_on_bdy(:,2) ./ Norm_on_bdy; % normalize
% update boundary points in global array (all other nodal values are unaffected)
Geo_q(dGamma_DoF_Indices,:) = Geo_q_on_bdy;
```
Thus, `Geo_q` is now a piecewise quadratic approximation of the unit disk.  Next,
```
% apply parameterization
Geo_Points = Psi(Geo_q(:,1), Geo_q(:,2));
% END: define the higher order surface
```

## Define Finite Element Spaces

```
% define FE spaces
V_Space    = FiniteElementSpace('V_h', P2_RefElem, Mesh, 'Gamma');
V_Space    = V_Space.Set_DoFmap(Mesh,uint32(G_DoFmap)); % G_Space has the same DoFmap

P1_RefElem = ReferenceFiniteElement(lagrange_deg1_dim1());
W_Space    = FiniteElementSpace('W_h', P1_RefElem, Mesh, 'dGamma');
W_DoFmap   = DEMO_LB_mex_W_Space_DoF_Allocator(uint32(Mesh.Get_Global_Subdomain('dGamma')));
W_Space    = W_Space.Set_DoFmap(Mesh,uint32(W_DoFmap));
```
we use the same DoFmap for `V_Space`, because `G_Space` and `V_Space` have the same reference finite element.  Note that we have a separate DoF allocator for `W_Space`.

## Assemble the Matrices

Because `Gamma` and `dGamma` are interacting sub-domains, we must compute the (topological) embedding data for the mesh:
```
Domain_Names = {'Gamma'; 'dGamma'};
Gamma_Embed = Mesh.Generate_Subdomain_Embedding_Data(Domain_Names);
```

Then (assuming we have mex'ed our bilinear form definition file), we assemble the matrices:
```
FEM = DEMO_mex_MatAssem_Laplace_Beltrami_Open_Surface([],Geo_Points,G_Space.DoFmap,[],...
                        Gamma_Embed,V_Space.DoFmap,W_Space.DoFmap);
```
Note that we feed the higher order mesh information `Geo_Points`, `G_Space.DoFmap` into the matrix assembler.

Extracting the matrices is easy with
```
% put FEM into a nice object to make accessing the matrices easier
LB_Mats = FEMatrixAccessor('Laplace-Beltrami',FEM);

% pull out the matrices
M = LB_Mats.Get_Matrix('M');
K = LB_Mats.Get_Matrix('K');
B = LB_Mats.Get_Matrix('B');
```

## Simple Error Checks

We numerically compute the length of `dGamma` using the `B` matrix:
```
dGamma_Subdomain = Mesh.Output_Subdomain_Mesh('dGamma');
disp('Check length of boundary curve:');
Len_mesh = sum(dGamma_Subdomain.Volume());
Len_by_matrix = sum(B(:));
```
and we compare to the "true" value computed analytically:
```
disp('Length by mesh: error is (small):');
Len_mesh - 7.640395578055425e+00
disp('Length by B matrix: error is (small):');
Len_by_matrix - 7.640395578055425e+00
```

## Define Finite Element Functions

Now we interpolate a known exact solution, and the corresponding right-hand-side data, using the coordinates of the DoFs in our *iso-parametric* FE space.  To do this, we need both `G_Space` and `V_Space`.  In fact, we need to do interpolation within `G_Space`.  Hence, we need to set a directory to store the mex file for the interpolation code:
```
G_Space = G_Space.Set_mex_Dir(Main_Dir,'mex_Laplace_Beltrami_Open_Surface_G_h_Interpolation');
```
Here, we named the "internal" interpolation mex file of `G_Space`.

Now we can get the precise coordinates of the DoFs in `V_Space`.
```
V_Points = V_Space.Get_DoF_Coord(Mesh,G_Space,Geo_Points);
```
Note: `V_Points` and `Geo_Points` are actually the same because `V_Space` and `G_Space` have the same reference finite element (P2 Lagrange) and the same DoFmap.

With this, we can interpolate given functions in order to compute the Right-Hand-Side `RHS` data:
```
f_h = f_exact(V_Points(:,1),V_Points(:,2));
g_h = u_exact(V_Points(:,1),V_Points(:,2));
```

## Solve the Linear System

Next, form the big system and solve:
```
W_Num = W_Space.num_dof;
ZZ = sparse(W_Num,W_Num);
MAT = [K, B';
       B, ZZ];
RHS = [M*f_h; B*g_h];
% use backslash to solve
Soln = MAT \ RHS;
```
Then parse the solution vector:
```
V_Num = V_Space.num_dof;
u_h = Soln(1:V_Num,1);
lambda_h = Soln(V_Num+1:end,1);
```

## Interpolate Normal and Tangent Vectors

```
% define exact normal and tangent vectors on dGamma
NV_func = @(u,v) [-u ./ sqrt(u.^2 + v.^2 + 1), v ./ sqrt(u.^2 + v.^2 + 1), 1 ./ sqrt(u.^2 + v.^2 + 1)];
TV_func = @(u,v) [-v ./ sqrt(u.^2 + v.^2 + 4*u.^2.*v.^2), u ./ sqrt(u.^2 + v.^2 + 4*u.^2.*v.^2), -2*u.*v ./ sqrt(u.^2 + v.^2 + 4*u.^2.*v.^2)];
% interpolate normal and tangent vectors using the G_Space nodal coordinates
NV_h = NV_func(Geo_Points(:,1),Geo_Points(:,2));
TV_h = TV_func(Geo_Points(:,1),Geo_Points(:,2));
```

## Compute Numerical Errors

Now (assuming we have mex'ed our error definition file), we assemble the errors:
```
FEM = DEMO_mex_Compute_Errors_Laplace_Beltrami_Open_Surface([],Geo_Points,G_Space.DoFmap,[],...
                        Gamma_Embed,G_Space.DoFmap,V_Space.DoFmap,W_Space.DoFmap,...
                        NV_h,TV_h,lambda_h,u_h);
```

Extracting the errors is easy with
```
% put FEM into a nice object to make accessing the matrices easier
Error_MATS = FEMatrixAccessor('Errors',FEM);

% pull out the (squared) errors, and take sqrt!
u_L2_Error = sqrt(Error_MATS.Get_Matrix('u_L2_Error_sq'));
lambda_L2_Error = sqrt(Error_MATS.Get_Matrix('lambda_L2_Error_sq'));
NV_Error = sqrt(Error_MATS.Get_Matrix('NV_Error_sq'));
TV_Error = sqrt(Error_MATS.Get_Matrix('TV_Error_sq'));
Neumann_L2_Error = sqrt(Error_MATS.Get_Matrix('Neumann_L2_Error_sq'));
```

And record each error in an array:
```
u_Error_Vec(kk) = u_L2_Error;
lambda_Error_Vec(kk) = lambda_L2_Error;
NV_Error_Vec(kk) = NV_Error;
TV_Error_Vec(kk) = TV_Error;
Neumann_Error_Vec(kk) = Neumann_L2_Error;
```

Finally, we close the *refinement* `for` loop:
```
end
```

# Plotting

We will plot the solution data corresponding to the finest mesh that was used.

## Process P2 Solution into P1 Solution

We cannot plot `u_h` directly because it is a function in `V_Space`.  We need to interpolate it at the vertices of the piecewise linear mesh.  Here is one way to do it that avoids writing an interpolation code:
```
% find the mapping from mesh points to P2 nodes
BB = [-1.001, 1.001, -1.001, 1.001, -1.001, 1.001];
OT = mexOctree(V_Points,BB);

% find the mapping via closest point
[P1_to_P2, OT_dist] = OT.kNN_Search(Mesh.Points,1);

% interpolate solution back onto piecewise linear mesh for plotting
u_h_P1 = u_h(P1_to_P2,1);

delete(OT); % clear the octree object!
```
Note that `mexOctree` is included in FELICITY.

## Plot `u_h_P1`

```
% plot u solution
FH_u_h = figure('Renderer','painters','PaperPositionMode','auto','Color','w');
trisurf(Mesh.ConnectivityList,Mesh.Points(:,1),Mesh.Points(:,2),Mesh.Points(:,3),u_h_P1);
shading interp;
colormap(jet);
lightangle(80,-20);
lightangle(-40,40);
lighting phong;
colorbar;
caxis([0 max(u_h)]);
title('Color plot of u_h on \Gamma');
AX = [-1 1 -1 1 -0.5 0.5];
axis(AX);
axis equal;
axis(AX);
view([-35,30]);
```

Here is a plot of `u_h_P1`:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Soln_u_h.jpg" width="600" alt="Laplace-Beltrami Solution u_h">

## Plot `lambda_h`

```
% Plot the \lambda solution:
Bdy_X = W_Space.Get_DoF_Coord(Mesh,G_Space,Geo_Points);
FH_lambda_h = figure('Renderer','painters','PaperPositionMode','auto','Color','w');
p2 = edgemesh(W_Space.DoFmap, Bdy_X, lambda_h);
set(p2,'LineWidth',2.0);
shading interp;
title('Color plot of \lambda_h on \partial \Gamma');
axis equal;
colorbar;
view([-35,30]);
```
Note that `edgemesh` is part of FELICITY.

Here is a plot of `u_h_P1`:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Soln_lambda_h.jpg" width="600" alt="Laplace-Beltrami Bdy Solution lambda_h">

## Plot The Error Decay

Now plot the error decay in order to view the convergence rate:
```
Order_2_Line = 4 * 4.^(-Refine_Vec);
Order_3_Line = 1.0 * 8.^(-Refine_Vec);
FH_error = figure('Renderer','painters','PaperPositionMode','auto','Color','w');
semilogy(Refine_Vec,u_Error_Vec,'k-*',Refine_Vec,lambda_Error_Vec,'r-.',...
         Refine_Vec,NV_Error_Vec,'b-d',Refine_Vec,TV_Error_Vec,'g-s',...
         Refine_Vec,Neumann_Error_Vec,'r--o',...
         Refine_Vec,Order_2_Line,'m-d',...
         Refine_Vec,Order_3_Line,'k-d');
xlabel('Refinement Level');
ylabel('L^2 errors');
title('Numerical Convergence Rates');
axis([2 7 1e-8 1e1]);
legend({'$\| u - u_h \|_{L^2(\Gamma)}$', '$\| \lambda - \lambda_h \|_{L^2(\Gamma)}$', ...
        '$\| \mathbf{\nu} - \boldmath{\nu}_h \|_{L^2(\Gamma)}$', '$\| \mathbf{\tau} - \mathbf{\tau}_h \|_{L^2(\Gamma)}$',...
        '$\| \lambda_h - (-\mathbf{\xi}_h \cdot \nabla_{\Gamma} u_h) \|_{L^2(\Gamma)}$',...
        '$O(h^2)$ line', '$O(h^3)$ line'},'Interpreter','latex','FontSize',12,'Location','best');
grid;
```

Here is a plot of the convergence rate:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Beltrami_Open_Surface_Conv_Rate.jpg" width="600" alt="Laplace-Beltrami Conv. Rate">

# Conclusion

The method clearly converges from the numerical tests.
