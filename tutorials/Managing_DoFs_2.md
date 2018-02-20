Managing Degrees-of-Freedom (DoFs): Part 2
==========================================

This tutorial describes how to correctly define the DoFmap for finite element spaces that are defined on *sub-domains* of a global mesh.

# Defining Finite Element Spaces on Sub-domains

Degree-of-Freedom (DoF) Allocation on the diagonal of a square:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Tutorial_Managing_DoFs_Ex_P1_On_1D_Subdomain.jpg|width=520|alt=DoF Allocation]]

## Introduction

Finite element spaces can be *defined over sub-domains* of a global mesh.  For example, you may have a piecewise linear space on a domain that is a polygonal curve embedded in a triangular mesh.  Allocating DoFs in this case is no more difficult than what was done in [Managing DoFs: Part 1](../wiki/Managing_DoFs_1).  But there may be some confusion between global mesh vertex indices and FE space DoF indices.  This tutorial explains the difference.

## Scalar-Valued Finite Element Space

The above figure is a simple mesh of a square domain with a piecewise linear finite element space defined over the diagonal.  The nodal variables are denoted by diamonds.  There are two local basis functions on each edge of the diagonal, so there are two indices associated with each edge:
```matlab
E_1 : 1 2
E_2 : 2 3
```
This is the DoFmap for the finite element space defined on the diagonal sub-domain of the mesh in the above figure.  In this case, the total number of global DoFs is three.

See Section 6.4.3 of the PDF manual for more info.

## Allocating DoFs On An Embedded Sub-domain

Let’s do an example. First create the mesh in the figure above by typing the following at the MATLAB prompt:
```matlab
% define a simple square mesh
Vtx = [0 0; 1 0; 1 1; 0 1; 0.5 0.5];
Tri = [1 2 5; 2 3 5; 3 4 5; 4 1 5];
% create a mesh object
Mesh = MeshTriangle(Tri,Vtx,'Square');
```

Now define the "diagonal" sub-domain of the square:
```matlab
% define a subdomain to be one of the diagonals of the square
Mesh = Mesh.Append_Subdomain('1D','Diag',[1 5; 5 3]);
```

Next, create a FiniteElementSpace object for the finite element space over the subdomain `'Diag'`:
```matlab
% declare reference element
P1_Elem = lagrange_deg1_dim1();
RE = ReferenceFiniteElement(P1_Elem);
% declare a finite element space that is defined on ’Diag’
FES = FiniteElementSpace('P1',RE,Mesh,'Diag');
clear RE;
```
Note that the reference element is 1-D, because Diag has topological dimension 1.

Now, set the DoFmap for this space. Since the domain consists of only two edges and the space is piecewise linear, we type the following commands:
```matlab
% set DoFmap
DFM = [1 2;
       2 3];
FES = FES.Set_DoFmap(Mesh,DFM);
clear DFM;
```
where edges `E_1` and `E_2` share DoF #2 at their common vertex.

Note: if the mesh has lots of elements, then you may want to automatically generate the code for creating the DoFmap (see [Allocating Degrees-of-Freedom (DoFs)](../wiki/Allocate_DoFs_1)).

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
text(XC(ii,1)+0.03,XC(ii,2)+0.0,DoFstr,'FontSize',16);
end
text(0.5,0.2,'T1');
text(0.8,0.5,'T2');
text(0.5,0.8,'T3');
text(0.2,0.5,'T4');
text(0.27,0.25,'E1');
text(0.77,0.75,'E2');
Mesh.Plot_Subdomain('Diag');
hold off;
axis equal;
axis(AX);
axis off;
title('P_1 DoFs of Finite Element Space (Numbered Diamonds)');
```

## Remarks On DoF Numbering

We emphasize that the DoF indices in the finite element space have *no relation* to the numbering of the vertices in the mesh. In this example, the center mesh vertex of the square has index #5, but the finite element space DoF index at the center vertex is #2.

When defining any DoFmap, there must not be any "gaps" in the DoF numbering in the DoFmap.  For example, we could not use this:
```matlab
% set DoFmap
DFM = [1 5;
       5 3];
FES = FES.Set_DoFmap(Mesh,DFM);
```
where `DFM` is simply the connectivity of Diag referenced to the global mesh vertices.  If we use this, then FELICITY will give an error message when calling `Set_DoFmap`.

In conclusion, if `N` is the number of distinct basis functions in the finite element space (i.e. the _dimension_ of the finite element space), then the DoFmap must contain *all indices* from 1 to N.