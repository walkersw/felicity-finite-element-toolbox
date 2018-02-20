Find Closest Points on a Surface Mesh
=====================================

This tutorial on how to find the closest points on a surface triangulation to a given set of arbitrary points in 3-D.

# Introduction

Finding the closest point to a manifold is a common task.  It is needed for interpolating data on the surface, or for computing distance functions.

FELICITY provides code generation to do this.  In fact, the underlying mesh geometry may be higher order (e.g. quadratic curved triangles).  The generated code is specific to the type of mesh geometry.

# Generating The Code

## Input File

First, generate code to search surface meshes in 3-D.  In the MATLAB editor, create the following m-function and name it `Point_Search_Surface.m`:

```matlab
function DOM = Point_Search_Surface()

% define domain (surface embedded in 3-D)
Surface = Domain('triangle',3);

% define geometry representation of global domain:  Domain, reference element
% Here we assume the default reference element:
% Lagrange piecewise linear continuous,
% in other words, the triangulation is flat (not curved).
G1 = GeoElement(Surface);

% define a set of domains to point search in
DOM = PointSearches(G1);

% collect together all of the domains to be searched
DOM = DOM.Append_Domain(Surface);

end
```

Don't worry about what it means just yet (the PDF manual explains more; see Chapter XX (to be added)). 

## Compile It

Put the file `Point_Search_Surface.m` into a directory that is *in your MATLAB path*.  Now compile it by typing the following command at the MATLAB prompt and press "ENTER":

```matlab
Convert_PtSearch_Definition_to_MEX(@Point_Search_Surface,{},'mex_Point_Search');
```

Here we named the executable `mex_Point_Search` (see next section).

# Computing Closest Points To A Sphere

A sphere may seem too easy of an example (since the closest points can be computed _exactly_ by hand).  However, it gives us a way to check the accuracy of the code.  And, the tutorial below works *exactly the same* for any other surface mesh.

## Initialization

Either in an m-script or at the MATLAB prompt, create a unit sphere triangle mesh:
```matlab
Refine_Level = 4;
Center = [0, 0, 0];
Radius = 1;
[TRI, VTX] = triangle_mesh_of_sphere(Center,Radius,Refine_Level);
Mesh = MeshTriangle(TRI, VTX, 'Surface');
```
The function `triangle_mesh_of_sphere` is a FELICITY routine.

## Create Points To Search

Define points in the global space:
```matlab
% define a regular grid of points
s_vec = (-1:0.6:1)'; % make sure no point is *exactly* at the origin
[XX,YY,ZZ] = meshgrid(s_vec,s_vec,s_vec);
GX = [XX(:), YY(:), ZZ(:)];
```

In order to search for the closest point on the surface triangulation, we need an *initial guess* for the particular triangle (i.e. cell) that contains the closest point.

The closer the initial guess, the faster the algorithm runs.  *Warning:* if the surface is very contorted and non-convex, and you give a bad initial guess, then the closest point computed by the MEX file may not be the true closest point.  The search algorithm may get trapped in a "local minimum."

For this example, we take the initial guess (for every point to be searched) to be triangle (cell) #1:
```matlab
Cell_Indices = uint32(ones(size(GX,1),1));
```

We also need a neighbor data structure for the mesh (to facilitate traversing the mesh).  This is done by the following command:
```matlab
Surface_Neighbors = uint32(Mesh.neighbors);
```

Next, put all this together into a cell array:
```matlab
Given_Points = {Cell_Indices, GX, Surface_Neighbors};
```
Note: you can also set `Cell_Indices` to an empty matrix.  In this case, the code assumes the initial guess to be cell #1 (for every point).

## Run The Search

Now we run the MEX file we generated before:
```matlab
SEARCH = mex_Point_Search(Mesh.Points,uint32(Mesh.ConnectivityList),[],[],Given_Points);
```

The output `SEARCH` is a struct of the form:
```matlab
SEARCH.DATA = cell array
SEARCH.Name = 'Surface'
```
where `SEARCH.DATA` is a 1x2 cell array.  The contents are:
```matlab
CI = double(SEARCH.DATA{1});
Local_Ref_Coord = SEARCH.DATA{2};
```
where `CI` are the cell (triangle) indices (found by the algorithm) that contain the closest points.

In addition, the coordinates of the closest point *with respect to the enclosing cell* are given in `Local_Ref_Coord`.  In other words, these coordinates are relative to the standard reference triangle.  This is very useful if you want to *interpolate* finite element functions defined on the surface mesh.

## Plot Results

Of course, we can convert these local coordinates to coordinates in the global space:
```matlab
% get global coordinates of the "found" closest points on sphere
XC = Mesh.referenceToCartesian(CI,Local_Ref_Coord);
```

And we can now plot it:
```matlab
FH = Mesh.Plot;
set(FH,'facealpha',0.4);
hold on;
plot3(GX(:,1),GX(:,2),GX(:,3),'k*','MarkerSize',10);
plot3(XC(:,1),XC(:,2),XC(:,3),'r.','MarkerSize',20);
plot3([XC(:,1), GX(:,1)]',[XC(:,2), GX(:,2)]',[XC(:,3), GX(:,3)]',...
      'b-','LineWidth',1.6);
hold off;
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
axis equal;
title('Closest Points On Triangulation Of Sphere');
```

This tutorial can be found in ".\Demo\Closest_Point_Sphere" sub-directory of FELICITY.  Also, see the PDF manual for another tutorial on finding points in planar triangulations.