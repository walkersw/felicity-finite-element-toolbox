Managing Degrees-of-Freedom (DoFs): Part 1
==========================================

This is a simple tutorial on how to manage and sift through DoFs using the FiniteElementSpace class in FELICITY.

# Managing Degrees-of-Freedom

Consider a Degree-of-Freedom (DoF) allocation on a simple two triangle mesh:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Tutorial_Managing_DoFs_Ex_P2.jpg|width=320|alt=DoF Allocation]]

## Scalar-Valued Finite Element Space

The above figure is a simple mesh of a square domain with a piecewise quadratic finite element space defined over it.  The nodal variables are denoted by diamonds.  There are six local basis functions on each triangle, so there are six indices associated with each triangle:

```matlab
T_1 : 1 2 3 8 6 5
T_2 : 1 3 4 9 7 6
```

This is the so-called Degree-of-Freedom map (*DoFmap*) for the finite element space defined on the mesh in the above figure. In this case, the total number of global DoFs is nine.

DoFmaps are stored as matrices in MATLAB. Each row of the DoFmap corresponds to a mesh element (e.g. triangle), and each column corresponds to the local DoFs for that element.  For example, column #2 of the DoFmap gives the global DoF indices for local basis function #2 of all the elements.

See Section 6.4 of the PDF manual for more info.

## Tensor-Valued Finite Element Space

For a tensor-valued space (such as when you are approximating a vector field), the DoFmap is a little more complicated.  Suppose the functions in our finite element space are vector-valued with two components.  Then the DoFmap decomposes into two parts:

```matlab
% DoFmap for component #1:
T_1 : 1 2 3 8 6 5
T_2 : 1 3 4 9 7 6

% DoFmap for component #2:
T_1 : 10 11 12 17 15 14
T_2 : 10 12 13 18 16 15
```
i.e. the DoFmap for component #2 is obtained by adding the number of DoFs in the scalar finite element space to each index of the DoFmap of component #1.

Thus, for storage purposes, we only need the DoFmap for the first component. The DoFmap for the second component can always be built on-the-fly from the first when needed.

See Section 6.4 of the PDF manual for more info.

# The FiniteElementSpace Class

FELICITY has the class FiniteElementSpace that allows for access to and manipulation of the DoFmap.

## Basic DoF Routines

Let’s do an example. First create the mesh in the figure above by typing the following at the MATLAB prompt:
```matlab
% define a simple square mesh
Vtx = [0 0; 1 0; 1 1; 0 1];
Tri = [1 2 3; 1 3 4];
% create a mesh object
Mesh = MeshTriangle(Tri,Vtx,'Square');
```

Now define several sub-domains of the square (see the mesh chapter 2 in the PDF manual):
```matlab
% define the boundary of the square
Bdy_Edge = Mesh.freeBoundary;
Mesh = Mesh.Append_Subdomain('1D','Boundary',Bdy_Edge);
% define some more subdomains
Mesh = Mesh.Append_Subdomain('0D','Bottom Left Corner',[1]);
Mesh = Mesh.Append_Subdomain('0D','Bottom Right Corner',[2]);
Mesh = Mesh.Append_Subdomain('0D','Top Right Corner',[3]);
Mesh = Mesh.Append_Subdomain('0D','Top Left Corner',[4]);
Mesh = Mesh.Append_Subdomain('1D','Bottom Side',[1 2]);
Mesh = Mesh.Append_Subdomain('1D','Top Side',[3 4]);
Mesh = Mesh.Append_Subdomain('1D','Diag',[1 3]);
Mesh = Mesh.Append_Subdomain('2D','Bottom Tri',[1]);
```

Next, load up the degree 2, Lagrange element of dimension 2 (i.e. defined on a triangle), and create a ReferenceFiniteElement object:
```matlab
% declare reference element
P2_Elem = lagrange_deg2_dim2();
RE = ReferenceFiniteElement(P2_Elem);
```

Now create the FiniteElementSpace object for the FE space, which is a cartesian (tensor) product of two scalar valued FE spaces:
```matlab
% declare a finite element space
FES = FiniteElementSpace('P2',RE,Mesh,[],2); % 2 components
clear RE; % not needed anymore
```

The arguments of FiniteElementSpace (in order) are:

* A name for the finite element space, e.g. `'P2'`.
* The reference finite element for the space.
* The mesh data for the finite element space.
* The name of the sub-domain that the finite element space is defined on. In the case above, the space is defined on the entire mesh (no sub-domain), so just pass the empty matrix [].
* The number of components of the space, if the space is a cartesian (tensor) product of scalar valued spaces.  If this argument is omitted, then the default value is 1 component.

We must still set the DoFmap for this space. So we store it by the following commands:
```matlab
% set DoFmap
DFM = [1 2 3 8 6 5;
       1 3 4 9 7 6];
FES = FES.Set_DoFmap(Mesh,DFM);
clear DFM; % not needed anymore
```
Note: if the mesh has lots of elements, then you may want to automatically generate the code for creating the DoFmap (see [Allocating Degrees-of-Freedom (DoFs)](../wiki/Allocate_DoFs_1)).

Now lets look at the object. We can get a summary about the finite element space by typing:
```matlab
FES
```
We can look at the DoFmap:
```matlab
FES.DoFmap
```
We can get a list of the DoF indices for the first component:
```matlab
DoF_Indices = FES.Get_DoFs
```

And we can get the global coordinates of the DoFs (for the first component):
```matlab
XC = FES.Get_DoF_Coord(Mesh)
```
Note that the mesh geometry is needed to get the coordinates.

We can plot the mesh and DoFs from the figure by the following commands:
```matlab
figure;
Mesh.Plot;
AX = [-0.05 1.05 -0.05 1.05];
axis(AX);
hold on;
for ii = 1:size(XC,1)
plot(XC(ii,1),XC(ii,2),'kd','LineWidth',1.7);
DoFstr = [num2str(ii)];
text(XC(ii,1)+0.03,XC(ii,2)+0.02,DoFstr,'FontSize',16);
end
text(0.75,0.25,'T1');
text(0.25,0.75,'T2');
hold off;
axis equal;
axis(AX);
axis off;
title('P_2 Degrees-of-Freedom (Numbered Diamonds)');
```

The DoFs for the second component of the FE space can be retrieved by:
```matlab
% get the 2nd component DoF indices
FES.Get_DoFs(2)
```
We can also get a matrix of DoF indices, where the k-th column is a list of the DoF indices for the k-th component:
```matlab
FES.Get_DoFs('all')
```
Note: each row of this matrix corresponds to one of the diamonds in the figure, i.e. each row stores the first and second component DoF indices, both of which have the same global coordinates.

## DoFs on Sub-domains

When sub-domains are present, you can retrieve the subset of DoFs that lie on a given sub-domain. For example, we can get the Boundary DoFs by:
```matlab
% get list of the (1st component) DoFs on the Boundary
Bdy_DoFs = FES.Get_DoFs_On_Subdomain(Mesh,'Boundary')
```

You can also get the DoFs for all tensor components on a given sub-domain by passing an additional argument `'all'`, e.g.
```matlab
% get list of all the (tensor component) DoFs on the Boundary
FES.Get_DoFs_On_Subdomain(Mesh,'Boundary','all')
```

Having defined `XC` in the previous section, we can get the coordinates of the DoFs on a sub-domain. For example, the coordinates of the DoFs on `'Boundary'` are obtained by:
```matlab
% get the coordinates of the Boundary DoFs
XC(Bdy_DoFs,:)
```

We can also store in the FES object which sub-domains have their DoFs fixed (i.e. for Dirichlet conditions). For example, we can fix the `'Bottom Side'` and `'Top Side'` sub-domains (defined earlier):
```matlab
FES = FES.Append_Fixed_Subdomain(Mesh,'Bottom Side');
FES = FES.Append_Fixed_Subdomain(Mesh,'Top Side');
```
Note: you must pass the mesh because the FiniteElementSpace object does some internal consistency checks to make sure the mesh actually contains the sub-domain.

You can now list the DoFs that are fixed by:
```matlab
FES.Get_Fixed_DoFs(Mesh)
FES.Get_Fixed_DoFs(Mesh,2) % if you want the fixed DoFs for component 2
```

Similarly, you can get the DoFs that are free by:
```matlab
FES.Get_Free_DoFs(Mesh)
```
If you want to clear the list of fixed sub-domains, just type
```matlab
FES = FES.Set_Fixed_Subdomains(Mesh,{});
```

The argument { } is an empty cell array. Alternatively, you could pass a cell array of strings, where each string is a name of a sub-domain. This is another way to set the fixed sub-domains, instead of using Append_Fixed_Subdomain.