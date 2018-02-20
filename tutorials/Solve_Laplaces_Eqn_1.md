Solving Laplace's Equation in 1-D
=================================

# Introduction

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Laplace_Var_Form.jpg|width=650|alt=Variational Form of Laplace's Eqn]]

To solve this numerically, we need to create the discrete operators (i.e. matrices). The philosophy behind FELICITY is to first define the bilinear forms (variational representation of the differential operators) via a user specified input file. Then FELICITY *automatically generates* a special purpose C++/MEX file and compiles it. This becomes an executable file that you can call directly from MATLAB. 

# Input File

In the MATLAB editor, create the following m-function and name it `Example_1D.m`:

```matlab
function MATS = Example_1D()

% define domain (1-D)
Omega = Domain('interval');

% define piecewise linear finite element space
Scalar_P1 = Element(Omega,lagrange_deg1_dim1,1);

v = Test(Scalar_P1); % test function
u = Trial(Scalar_P1); % trial function

Stiff_Matrix = Bilinear(Scalar_P1,Scalar_P1); % stiffness matrix
I1 = Integral(Omega,v.grad' * u.grad);
Stiff_Matrix = Stiff_Matrix.Add_Integral(I1);

RHS = Linear(Scalar_P1);
RHS = RHS + Integral(Omega,-1.0 * v.val);

Quadrature_Order = 3; % order of accuracy for the quad rule
G1 = GeoElement(Omega); % define geometric representation of the domain

MATS = Matrices(Quadrature_Order,G1); % collect all matrices together
MATS = MATS.Append_Matrix(Stiff_Matrix);
MATS = MATS.Append_Matrix(RHS);

end 
```

Don't worry about what it means just yet (the PDF manual explains more). 

# Compile It!

Put the file `Example_1D.m` into a directory that is *in your MATLAB path*.  Now compile it by typing the following command at the MATLAB prompt and press "ENTER":

```matlab
Convert_Form_Definition_to_MEX(@Example_1D,{},'Assemble_1D');
```

Here we named the executable `Assemble_1D` (see next section).

# Run It!

We will solve the 1-D Laplace problem. First, type the following commands at the MATLAB prompt:

```matlab
>> Num_Vertices = 101;
>> Indices = (1:1:Num_Vertices)';
>> Omega_Mesh = uint32([Indices(1:end-1,1), Indices(2:end,1)]);
>> Omega_Vertices = linspace(0,1,Num_Vertices)';
```

This defines a mesh that partitions [0, 1] into 100 sub-intervals.

Next, we need the local-to-global Degree-of-Freedom map (DoF map) for the finite element space. Since we are using piecewise linear continuous basis functions, we can simply type:

```matlab
>> P1_DoFmap = Omega_Mesh;
```

Now assemble the matrix representations of the bilinear (and linear) forms by executing:

```matlab
>> FEM = Assemble_1D([],Omega_Vertices,Omega_Mesh,[],[],P1_DoFmap);
```

`FEM` is a MATLAB struct that contains the sparse matrices in alphabetical order based on the names used in the `Example_1D.m` file. Finally, we obtain the solution of the PDE by the following commands:

```matlab
>> A = FEM(2).MAT; % stiffness matrix

>> Soln = zeros(size(A,1),1); % init

>> alpha = 0;
>> RHS = FEM(1).MAT;
>> RHS(end) = RHS(end) + alpha; % add in the Neumann condition at s = 1

>> % impose Dirichlet condition at s=0
>> Soln(2:end) = A(2:end,2:end) \ RHS(2:end);
```

You can plot the solution with the following commands:

```matlab
>> figure;
>> h1 = plot(Omega_Vertices,Soln,'b-','LineWidth',2.0);
>> title('Solution of PDE: -d^2/ds^2 u = -1, u(0) = 0, d/ds u(1) = 0','FontSize',14);
>> xlabel('x (domain = [0, 1])','FontSize',14);
>> ylabel('u (solution value)','FontSize',14);
>> set(gca,'FontSize',14);
>> axis equal;
```

For more information, see the *Assembling Matrices* chapter in the PDF manual. 