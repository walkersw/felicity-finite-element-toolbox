Interpolating Finite Element Data
=================================

This is a tutorial on generating code to perform interpolation of finite element data.

# Introduction

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Tutorial_Interpolation_Ex_P2.jpg|width=800|alt=Interpolating P2 Finite Element Data]]

FELICITY can *automatically generate* a special purpose C++/MEX file and compile it for most interpolation purposes. This becomes an executable file that you can call directly from MATLAB. 

# Input File

In the MATLAB editor, create the following m-function and name it `Interpolate_Grad_P_X_2D.m`: 

```matlab
function INTERP = Interpolate_Grad_P_X_2D()

% define domain
Omega = Domain('triangle');

% define finite element spaces
Scalar_P2  = Element(Omega,lagrange_deg2_dim2,1);

% define functions on FE spaces
p = Coef(Scalar_P2);

% define geometric function on 'Omega' domain
gf = GeoFunc(Omega);

% define expressions to interpolate
I_grad_p_X = Interpolate(Omega,p.grad' * gf.X);

% define geometry representation - Domain, reference element
G1 = GeoElement(Omega);

% define a set of interpolations to perform
INTERP = Interpolations(G1);

% collect all of the interpolations together
INTERP = INTERP.Append_Interpolation(I_grad_p_X);

end
```

Don't worry about what it means just yet (the PDF manual explains more; see Chapter 7). 

# Compile It!

Put the file `Interpolate_Grad_P_X_2D.m` into a directory that is *in your MATLAB path*.  Now compile it by typing the following command at the MATLAB prompt and press "ENTER":

```matlab
Convert_Interp_Definition_to_MEX(@Interpolate_Grad_P_X_2D,{},'Interp_2D');
```

Here we named the executable `Interp_2D` (see next section).

# Run It!

First, create the domain (triangulation) for the unit square. At the MATLAB prompt (or in a script) type the following commands:
```matlab
% create mesh (Omega is the unit square)
Vtx = [0 0; 1 0; 1 1; 0 1];
Tri = [1 2 3; 1 3 4];
Mesh = MeshTriangle(Tri, Vtx, 'Omega');
```
i.e. the mesh consists of two triangles. Next, define the piecewise quadratic finite element space:
```matlab
% define FE space
P2_RefElem = ReferenceFiniteElement(lagrange_deg2_dim2());
P2_Lagrange_Space = FiniteElementSpace('Scalar_P2', P2_RefElem, Mesh, 'Omega');
P2_DoFmap = uint32(Setup_Lagrange_P2_DoFmap(Mesh.ConnectivityList,[]));
P2_Lagrange_Space = P2_Lagrange_Space.Set_DoFmap(Mesh,P2_DoFmap);
P2_X = P2_Lagrange_Space.Get_DoF_Coord(Mesh); % get coordinates of P2 DoFs
```
See [Managing DoFs: Part 1](../wiki/Managing_DoFs_1) for more info on how the FiniteElementSpace class works.

Set the coefficient function `p(x, y) = sin(x + y)`:
```matlab
% define coefficient function
p_func = @(x,y) sin(x + y);
px_func = @(x,y) cos(x + y);
py_func = @(x,y) cos(x + y);
p_val = p_func(P2_X(:,1),P2_X(:,2)); % coefficient values
```

Next, define an analytic function for the expression ∇p · x:
```matlab
% exact interpolant function
I_grad_p_X = @(x,y) px_func(x,y) .* x + py_func(x,y) .* y;
```

## Interpolate At Two Points

Now, define the coordinates of the interpolation points. In this example, we will only have two interpolation points:
```matlab
% define interpolation points

% specify the triangle indices to interpolate within
Tri_Indices = [1;
               2];
% specify interpolation coordinates w.r.t. reference triangle
Ref_Coord = [(1/3), (1/3);
             (1/3), (1/3)];
% collect into a cell array
Omega_Interp_Data = {uint32(Tri_Indices), Ref_Coord};
```

The format requires the interpolation point data to be stored in a MATLAB cell array:
* The first element of Omega_Interp_Data is a MATLAB vector of unsigned integers, of length `R`, which are triangle indices. Note: in order to interpolate, you must specify the mesh elements to interpolate within.
* The second element of Omega_Interp_Data is a MATLAB matrix of size `R×T`; `R` is the number of interpolation points and `T` is the topological dimension of the mesh element. In this example, `T = 2`.

Now use the executable that was generated before. Run the following at the MATLAB prompt (or insert after the above code in a script):
```matlab
INTERP = Interp_2D(Mesh.Points,uint32(Mesh.ConnectivityList),[],[],Omega_Interp_Data,P2_DoFmap,p_val);
```
Recall that we named the MEX file `Interp_2D`. In general, if you try to run the MEX file with the incorrect number of arguments, it will give you an error message but it will also tell you what the inputs should be.

`INTERP` is a MATLAB struct that contains the interpolation data in alphabetical order based on the names used in the `Interpolate_Grad_P_X_2D.m` file:
```matlab
INTERP(1).Name : 'I_grad_p_X'
INTERP(1).DATA : cell array containing interpolation data for the expression
```

## Interpolate At Many Points

Interpolating two points is not so complicated. Lets create a grid of points on the unit square and interpolate at those:
```matlab
% now interpolate at a lot more points
x_vec = (0:0.1:1);
y_vec = (0:0.1:1);
[XX, YY] = meshgrid(x_vec,y_vec);
Interp_Pts_2 = [XX(:), YY(:)];
```

Remember that interpolating at a given point requires we know the mesh element that the point belongs to. Thus, we must identify which triangle the above grid points belong to:
```matlab
% find which triangle the points belong to
BOT_TRI = (Interp_Pts_2(:,1) >= Interp_Pts_2(:,2));
Tri_Indices_2 = zeros(size(Interp_Pts_2,1),1);
Tri_Indices_2(BOT_TRI) = 1; % bottom triangle cell index
Tri_Indices_2(~BOT_TRI) = 2; % top triangle cell index
```

Then, we must convert the global grid coordinates to local reference coordinates:
```matlab
Ref_Coord_2 = Mesh.cartesianToReference(Tri_Indices_2,Interp_Pts_2);
Omega_Interp_Data_2 = {uint32(Tri_Indices_2), Ref_Coord_2};
```

Execute the interpolation code again:
```matlab
% interpolate!
INTERP_2 = Interp_2D(Mesh.Points,uint32(Mesh.ConnectivityList),[],[],Omega_Interp_Data_2,P2_DoFmap,p_val);
Interp_Values = INTERP_2(1).DATA{1,1};
```

Now plot the grid data:
```matlab
figure;
surf(XX, YY, I_grad_p_X(XX,YY)); % plot "exact" surface
hold on;
% plot the FE interpolated data as black dots
plot3(Interp_Pts_2(:,1),Interp_Pts_2(:,2),Interp_Values,'k.','MarkerSize',18);
hold off;
axis equal;
AZ = 70;
EL = 10;
view(AZ,EL);
title('Surface is Exact. Points Are Interpolated From FE Approximation','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12);
set(gca,'FontSize',12);
```
Note that the black dots do not appear exactly on the surface. This makes sense because we are interpolating a finite element *approximation* of the surface.

For more information, see the *Interpolating Finite Element Data* chapter in the PDF manual.