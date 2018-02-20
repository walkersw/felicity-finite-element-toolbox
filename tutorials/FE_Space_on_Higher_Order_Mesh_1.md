A Finite Element Space on a Higher Order Mesh
=============================================

This tutorial describes how to do higher order meshes (e.g. piecewise quadratic triangles) instead of the standard piecewise linear case.  The details can be found in the FELICITY manual, the chapter on "Helper Routines For Finite Element Codes".

Note: this example can be run using the script `test_Curved_Domain_2D.m` located in
```
./FELICITY/Demo/Curved_Domain_2D/
```

# Domain Definition

We will define a piecewise quadratic triangulation of the following parameterized domain:
[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Curved_Domain_Param.jpg|width=320|alt=Curved Domain Param]]

A piecewise linear mesh approximation of the above domain is shown below:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Curved_Domain_Diagram.jpg" width="450" alt="Curved Domain Diagram">

# Example Finite Element Forms

We will define a continuous piecewise quadratic (Lagrange) finite element space on the piecewise quadratic mesh, a so-called iso-parametric FE space.  In addition, we will define the following bilinear and linear forms on the FE space:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Curved_Domain_Forms.jpg" width="300" alt="Curved Domain FE Forms">

where the last linear form computes the integral of the Neumann data of `v_h`.

# Abstract Definition File

The (abstract) m-file that defines all of this is given here:
```
function MATS = MatAssem_Curved_Domain_2D()

% define domain (in the x-y plane)
Omega = Domain('triangle');
Gamma = Domain('interval') < Omega; % subset

% define finite element spaces
V_h = Element(Omega, lagrange_deg2_dim2); % piecewise quadratic

% define functions in FE spaces
v_h = Test(V_h);
u_h = Trial(V_h);

% geometric function
gf = GeoFunc(Gamma);

% define (discrete) forms
M = Bilinear(V_h,V_h);
M = M + Integral(Omega, v_h.val * u_h.val );

M_Gamma = Bilinear(V_h,V_h);
M_Gamma = M_Gamma + Integral(Gamma, v_h.val * u_h.val );

B_neu = Linear(V_h);
B_neu = B_neu + Integral(Gamma, v_h.grad' * gf.N );

% set the minimum order of accuracy for the quad rule
Quadrature_Order = 10;
% define how \Omega is represented: with piecewise quadratics
G1 = GeoElement(Omega,lagrange_deg2_dim2);
% define a set of matrices
MATS = Matrices(Quadrature_Order,G1);

% collect all of the matrices together
MATS = MATS.Append_Matrix(B_neu);
MATS = MATS.Append_Matrix(M);
MATS = MATS.Append_Matrix(M_Gamma);

end
```
Note that the space `V_h` is piecewise quadratic, and the `GeoElement` is defined to be piecewise quadratic.

# Implementation

After compiling the abstract definition file (for example, see [Solve Laplace's Equation](../wiki/Solve_Laplaces_Eqn_1)), we now implement a script for assembling the previously defined forms (FE matrices) and computing some basic consistency checks.  This involves doing a convergence check, i.e. we will loop over several meshes of decreasing mesh size.

Note: this can be found in `Execute_Curved_Domain_2D.m`.

## Initialization

We start with some initialization:
```
Main_Dir = <choose a dir in your MATLAB path>;

Refine_Vec = [0 1 2 3 4 5];
Num_Refine = length(Refine_Vec);

% init two arrays, one for each mapping method
Neumann_Error_Vec = {zeros(Num_Refine,1); zeros(Num_Refine,1)};
```

Now define the domain parameterization, as well as parameterizations of the bottom and top boundary curves of the domain:
```
% define domain parameterization
Psi = @(u,v) [u, v + u.*(1-u).*(v-0.5)];
Gamma_y_bot = @(x) x.*(1-x).*(-0.5);
Gamma_y_top = @(x) 1 + x.*(1-x).*(1-0.5);
```

Next, start the big loop:
```
% mm==1 % square
% mm==2 % mapped square
for mm = 1:2
for kk = 1:Num_Refine
```
where we loop over two different methods for defining the piecewise quadratic domain; we also loop over successively finer meshes.

## Create The Base (Piecewise Linear) Mesh

Now create a mesh of the unit triangle `[0, 1] x [0, 1]`:
```
% BEGIN: define reference mesh
Refine_Level = Refine_Vec(kk);
Np = 2^Refine_Level + 1;
[Tri, Pts_u_v] = bcc_triangle_mesh(Np,Np);
```
and we either keep the unit square, or we map it (depending on the "method" chosen):
```
if (mm==1)
    Mesh_Pts = Pts_u_v;
else
    % apply parameterization
    Mesh_Pts = Psi(Pts_u_v(:,1), Pts_u_v(:,2));
end
Mesh = MeshTriangle(Tri, Mesh_Pts, 'Omega');
% END: define reference mesh
```
This gives us a base piecewise *linear* mesh.

Add the boundary `Gamma` as a sub-domain:
```
% add the open boundary
BDY  = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('1D','Gamma',BDY);
```

## Define a GeoElementSpace

Next, we create a special object for handling the higher order mesh geometry:
```
% create GeoElementSpace
P2_RefElem = ReferenceFiniteElement(lagrange_deg2_dim2());
G_Space  = GeoElementSpace('G_h',P2_RefElem,Mesh);
G_DoFmap = DEMO_Curved_Domain_mex_V_Space_DoF_Allocator(uint32(Mesh.ConnectivityList));
G_Space  = G_Space.Set_DoFmap(Mesh,uint32(G_DoFmap));
```
Note: at this point, we have already created a mex file to allocate DoFs for the piecewise quadratic space.

## Create the Piecewise Quadratic Mapping

With this, we get an initial piecewise quadratic FE function that represents the identity map for our mesh:
```
% BEGIN: define the higher order domain
Geo_Points_hat = G_Space.Get_Mapping_For_Piecewise_Linear_Mesh(Mesh);
% initially starts as a piecewise quadratic polynomial (P2 function) interpolating a
% piecewise linear (P1) function (representing the initial square mesh)
```
Since the base mesh is piecewise *linear*, this "quadratic" function is actually interpolating a piecewise linear function.

Now, depending on the "method" we chose, we create a new variable `Geo_Points` that represents the true piecewise quadratic mapping that we need to define our curved domain:
```
if (mm==1)
    % apply parameterization
    Geo_Points = Psi(Geo_Points_hat(:,1), Geo_Points_hat(:,2)); 
else
    DoFs_on_Gamma = G_Space.Get_DoFs_On_Subdomain(Mesh,'Gamma');
    GP_Gamma = Geo_Points_hat(DoFs_on_Gamma,:);
    Bot_Mask = GP_Gamma(:,2) < 0 + 1e-12; % <= zero
    Top_Mask = GP_Gamma(:,2) > 1 - 1e-12; % >= 1
    GP_Gamma(Bot_Mask,2) = Gamma_y_bot(GP_Gamma(Bot_Mask,1));
    GP_Gamma(Top_Mask,2) = Gamma_y_top(GP_Gamma(Top_Mask,1));
    Geo_Points = Geo_Points_hat;
    Geo_Points(DoFs_on_Gamma,:) = GP_Gamma;
    clear DoFs_on_Gamma GP_Gamma Bot_Mask Top_Mask;
end
% END: define the higher order domain
```
*Method 1* is the simplest, but it requires mapping from the unit square to the curved domain (so `Geo_Points` is not a small perturbation).  Furthermore, Method 1 has *all* of the elements being curved.

*Method 2* has the base mesh being a piecewise linear approximation of the curved domain; thus, the P2 function `Geo_Points` is a small perturbation of the piecewise linear approximation.  Indeed, only the triangles with an edge on the bottom or top boundary are curved; all other triangles are actually straight, even though a P2 mapping is used.

Depending on the situation, one method may be more useful than another.

## Define a Finite Element Space

```
V_Space = FiniteElementSpace('V_h', P2_RefElem, Mesh, 'Omega');
V_Space = V_Space.Set_DoFmap(Mesh,uint32(G_DoFmap)); % G_Space has the same DoFmap
```
we use the same DoFmap, because `G_Space` and `V_Space` have the same reference finite element.

## Assemble the Matrices

Because `Omega` and `Gamma` are interacting sub-domains, we must compute the (topological) embedding data for the mesh:
```
Domain_Names = {'Omega'; 'Gamma'};
Omega_Embed = Mesh.Generate_Subdomain_Embedding_Data(Domain_Names);
```

Then (assuming we have mex'ed our form definition file), we assemble the matrices:
```
FEM = DEMO_mex_MatAssem_Curved_Domain_2D([],Geo_Points,G_Space.DoFmap,[],...
                        Omega_Embed,V_Space.DoFmap);
```
Note that we feed the higher order mesh information `Geo_Points`, `G_Space.DoFmap` into the matrix assembler.

Extracting the matrices is easy with
```
% put FEM into a nice object to make accessing the matrices easier
my_Mats = FEMatrixAccessor('Curved Domain',FEM);

% pull out the matrices
M = my_Mats.Get_Matrix('M');
M_Gamma = my_Mats.Get_Matrix('M_Gamma');
B_neu = my_Mats.Get_Matrix('B_neu');
```

## Simple Error Checks

We numerically compute the area of the domain and the length of its perimeter:
```
% check domain measures
Area_by_matrix = sum(M(:));
Len_by_matrix = sum(M_Gamma(:));
```
and we compare to the "true" value computed analytically:
```
disp('Area by FE matrix: error is (small):');
Area_by_matrix - (7/6)
disp('Length by FE matrix: error is (small):');
Len_by_matrix - 4.080457638869102
```

## Define a Finite Element Function

Now we interpolate a known function, using the coordinates of the DoFs in our *iso-parametric* FE space.  To do this, we need both `G_Space` and `V_Space`.  In fact, we need to do interpolation within `G_Space`.  Hence, we need to set a directory to store the mex file for the interpolation code:
```
G_Space = G_Space.Set_mex_Dir(Main_Dir,'mex_Curved_Domain_G_h_Interpolation');
```
Here, we named the "internal" interpolation mex file of `G_Space`.

Now we can get the precise coordinates of the DoFs in `V_Space`.
```
V_Points = V_Space.Get_DoF_Coord(Mesh,G_Space,Geo_Points);
```
Note: `V_Points` and `Geo_Points` are actually the same because `V_Space` and `G_Space` have the same reference finite element (P2 Lagrange) and the same DoFmap.

With this, we can interpolate a given smooth function:
```
f_exact = @(x,y) sin(y);
f_h = f_exact(V_Points(:,1),V_Points(:,2));
```

## Consistency Check

We use `B_neu` to compute the integrated Neumann data for `f_h`, and we compare to the "exact" analytic value:
```
% compute Neumann data for known function
Neumann_data_exact = -5.277816286458646e-01;
Neumann_data_by_matrix = f_h' * B_neu;
Neumann_Error = abs(Neumann_data_exact - Neumann_data_by_matrix);
Neumann_Error_Vec{mm}(kk) = Neumann_Error;
```

Finally, we close the two `for` loops:
```
end

end
```

# Plotting

Plot the mesh of `Omega` and highlight boundary `Gamma`:
```
FH_mesh = figure('Renderer','painters','PaperPositionMode','auto','Color','w');
Mesh.Plot;
hold on;
Mesh.Plot_Subdomain('Gamma');
hold off;
title('Curved Domain \Omega (Piecewise Linear Approximation)');
%title('Om');
grid on;
axis([-0.2 1.2 -0.2 1.2]);
axis equal;
axis([-0.2 1.2 -0.2 1.2]);
```

Now plot the error decay in order to view the convergence rate:
```
Order_2_Line = 0.1 * 4.^(-Refine_Vec);
FH_error = figure('Renderer','painters','PaperPositionMode','auto','Color','w');
semilogy(Refine_Vec,Neumann_Error_Vec{1},'r-o',...
         Refine_Vec,Neumann_Error_Vec{2},'b-o',...
         Refine_Vec,Order_2_Line,'m-d');
xlabel('Refinement Level');
ylabel('Error');
title('Numerical Convergence Rate');
axis([0 5 1e-5 1e-1]);
legend({'$| Q - \int_{\Gamma} \nabla v_h \cdot \mathbf{\nu} |$ (method 1)',...
        '$| Q - \int_{\Gamma} \nabla v_h \cdot \mathbf{\nu} |$ (method 2)',...
        '$O(h^2)$ line'},'Interpreter','latex','FontSize',12,'Location','best');
grid;
```

Here is a plot of the convergence rate:

<img src="https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Curved_Domain_Conv_Rate.jpg" width="500" alt="Curved Domain Conv. Rate">

# Conclusion

One can use higher than degree 2 Lagrange elements for the geometry discretization.
