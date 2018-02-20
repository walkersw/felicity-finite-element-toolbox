Mesh Generation With TIGER: Part 2
==================================

This is a tutorial on generating unstructured meshes of polygons and surface meshes.

# Introduction

The TIGER algorithm can also mesh domains described by a closed polygonal curve or a water-tight polyhedral surface.

Note: the output mesh from the TIGER algorithm will not have the same vertices as the polygon or polyhedral surface used to describe the domain.  Also, you can only use bisection to compute the cut points.

# Example in 2-D

We shall use the class `Polygon_Mesh` located in
```matlab
./FELICITY/Static_Codes/Isosurface_Meshing/LevelSets_2D/@Polygon_Mesh
```
to create a 2-D mesh of the domain.

## Create Polygonal Mesh

First, sample a parametric curve:
```matlab
t = linspace(0,1,1001)';
x = (1 + 0.3*sin(5*2*pi*t)) .* cos(2*pi*t);
y = (1 + 0.3*sin(5*2*pi*t)) .* sin(2*pi*t);
```
Next, create a polygonal mesh:
```matlab
Vtx = [x,y];
Vtx = Vtx(1:end-1,:);
Ind = (1:1:size(Vtx,1))';
Edge = [Ind(1:end-1,1), Ind(2:end,1); Ind(end,1), Ind(1,1)];
```
Then, create a Level Set object:
```matlab
LS = Polygon_Mesh(Vtx,Edge);
```

## Mesh It!

Create the 2-D Mesher:
```matlab
Box_Dim = [-2, 2];
Num_BCC_Points = 100;
Use_Newton = false; % cannot use Newton Method here...
TOL = 1E-3;
MG = Mesher2Dmex(Box_Dim,Num_BCC_Points,Use_Newton,TOL);
```

Next, we need to pass a function handle for the level set interpolation to the mesh generator.  For this example, we do
```matlab
% setup up handle to interpolation routine
Interp_Handle = @(pt) LS.Interpolate(pt);
```
In general, the user may provide their own interpolation routine that samples the level set function.  The format of the function must be:
```matlab
[phi, grad_phi] = Interp_Func(point);
```
where
```matlab
% phi      = the level set function value
% grad_phi = the gradient of the level set function
% point    = array of point coordinate to evaluate at
```

Generate the mesh:
```matlab
MG = MG.Get_Cut_Info(Interp_Handle);
[TRI, VTX] = MG.run_mex(Interp_Handle);
```

## View It

Plot the result:
```matlab
p1 = trimesh(TRI,VTX(:,1),VTX(:,2),0*VTX(:,2));
set(p1,'EdgeColor','k');
axis equal;
grid off;
AX = 1.5*[-1 1 -1 1];
axis(AX);
```
It should look like this:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Star_Mesh_Polygon.jpg|width=500|alt=Star Mesh]]

# Example in 3-D

We shall use the class `Surface_Mesh` located in
```matlab
./FELICITY/Static_Codes/Isosurface_Meshing/LevelSets_3D/@Surface_Mesh
```
to create a 3-D mesh of the domain.

## Load Surface Mesh

For simplicity, we will use a surface mesh that is included with FELICITY:
```matlab
BT = load('Bumpy_Torus.mat','VTX','TRI');
LS = Surface_Mesh(BT.VTX,BT.TRI);
clear BT;
```

## Mesh It!

Create the 3-D Mesher:
```matlab
Box_Dim = [0, 1];
Num_BCC_Points = 41;
Use_Newton = false; % cannot use Newton Method here...
TOL = 1E-2;
MG = Mesher3Dmex(Box_Dim,Num_BCC_Points,Use_Newton,TOL);
```
Generate the mesh:
```matlab
% setup up handle to interpolation routine
Interp_Handle = @(pt) LS.Interpolate(pt);

MG = MG.Get_Cut_Info(Interp_Handle);
[TET, VTX] = MG.run_mex(Interp_Handle);
```

Extract the surface mesh of the tetrahedral mesh:
```matlab
% create mesh object
MT = MeshTetrahedron(TET, VTX, 'Bumpy Torus');
FACES = MT.freeBoundary;
```

## View It

Plot the result:
```matlab
p1 = trimesh(FACES,VTX(:,1),VTX(:,2),VTX(:,3));
set(p1,'EdgeColor','k');
axis equal;
grid off;
AX = [0 1 0 1 0 1];
axis(AX);
AZ = -80;
EL = 60;
view(AZ,EL);
grid on;
title('Surface of Bumpy Torus Mesh');
```
It should look like this:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Bumpy_Torus_Mesh_Surface.jpg|width=670|alt=Bumpy Torus]]

You can also plot it with some fancy lighting and shading:
```matlab
% make a nice plot of surface
figure;
p2 = patch('Vertices', MT.Points, 'Faces', FACES);
PURPLE = 2.5*[36 0 84]/255; % brighten it
set(p2,'FaceColor',PURPLE,'EdgeColor','none');
%set(p2,'FaceColor','magenta','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
AZ1 = -130;
EL1 = 30;
view(AZ1,EL1);
AZ2 = -80;
EL2 = 10;
camlight(AZ2,EL2);
lighting gouraud;
axis equal;
axis(AX);
```
It should look like this:

[[https://github.com/walkersw/felicity-finite-element-toolbox/blob/master/images/Bumpy_Torus_Mesh_Surface_With_Lighting.jpg|width=670|alt=Bumpy Torus With Lighting]]
