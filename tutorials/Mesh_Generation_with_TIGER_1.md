Mesh Generation With TIGER: Part 1
==================================

This is a tutorial on generating unstructured meshes of domains described by iso-surfaces.

# Introduction

FELICITY provides an algorithm for generating unstructured (triangle and tetrahedral) meshes with *guaranteed* angle bounds.  In 3-D, the dihedral angles are mathematically guaranteed to be between *8.54* degrees and *164.18* degrees.  Hence, meshes can be generated _robustly_.  And it is very *fast*!

*Caveat:* The method only creates meshes of _interior_ regions.  For example, suppose you have a square domain that contains a circular hole.  You can either mesh the part of the square that is outside the hole, or you can mesh the hole alone.  But you cannot mesh both simultaneously, such that the resulting mesh has edges that conform to the boundary of the hole.  If you try to do this, then the angle bounds are not guaranteed.  (Note: for 2-D meshes, it seems you *can* get a conforming mesh of both interior and exterior regions.)

Also note these meshes are quasi-uniform.

# Simple Example (2-D)

Mesh an ellipse.

## Initialization

After you have installed FELICITY, make sure you run `test_FELICITY.m` or `compile_static_codes.m`.  This will compile the C++ code that is the backbone of the mesh generator.

We start simply with an ellipse in 2-D.  Type the following at the MATLAB prompt (or put it in a script file):
```matlab
% create mesher object
Box_Dim = [0, 1];
Num_BCC_Points = 31;
Use_Newton = false;
TOL = 1E-2; % tolerance to use in computing the cut points
% BCC mesh of the unit square [0,1] x [0,1]
MG = Mesher2Dmex(Box_Dim,Num_BCC_Points,Use_Newton,TOL);
```

This creates a MATLAB object that runs the mesher and interfaces with a MEX/C++ code.  Here, we have created a background mesh with 31 points along the x and y dimensions for a unit square.  We have chosen not to use Newton's method for computing cut point intersections of the iso-surface with the mesh edges (default is simple bisection).

The tolerance dictates how well we approximate the cut points.  Here, it is set to `1E-2` and is relative to the mesh edge lengths. Thus it is a relative tolerance (or relative error).

## Define Iso-surface

Next, we define a level set function that represents an ellipse:
```matlab
LS = LS_Many_Ellipses();
LS.Param.rad_x = [0.4];
LS.Param.rad_y = [0.2];
LS.Param.cx    = [0.5];
LS.Param.cy    = [0.5];
LS.Param.sign  = [1]; % +1 or -1 only
```
which defines the two radii, the center position of the ellipse, and whether it is a positive or negative ellipse (inside or out).  Note: `LS_Many_Ellipses` is a simple MATLAB class (provided by FELICITY) that implements an interpolation routine that defines the iso-surface.  Look in the following directory for more info:
```matlab
./FELICITY/Static_Codes/Isosurface_Meshing/LevelSets_2D/@LS_Many_Ellipses
```

## Generate Mesh

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

Then, run the mesher with the following commands.
```matlab
MG = MG.Get_Cut_Info(Interp_Handle);
[TRI, VTX] = MG.run_mex(Interp_Handle);
```

This outputs the triangle connectivity data `TRI` (Mx3 matrix, where M is the number of triangles), and the vertex coordinates `VTX` (Nx2 matrix, where N is the number of vertices).

Now you can plot the mesh with the MATLAB command trimesh:
```matlab
trimesh(TRI,VTX(:,1),VTX(:,2));
axis equal;
grid on;
```
Note: you can also use the FELICITY class `MeshTriangle` to plot it and make other manipulations.  See the tutorial [Mesh Classes](../wiki/Tutorial_Meshes_1) for more info.

# Mesh With Holes (2-D)

Mesh a disk with ellipse shaped holes.

## Define Iso-surface

Use the same `LS_Many_Ellipses` class to define a domain with holes:
```matlab
LS = LS_Many_Ellipses();
LS.Param.rad_x = [0.45, 0.1, 0.2];
LS.Param.rad_y = [0.45, 0.2, 0.1];
LS.Param.cx    = [0.5, 0.35, 0.65];
LS.Param.cy    = [0.5, 0.35, 0.65];
LS.Param.sign  = [1, -1, -1];
```

All of the parameters are vectors of equal length that define three ellipses.  The first ellipse is actually a circle (equal x and y radii) and has +1 sign.  The next two ellipses have a negative sign which indicates to subtract them (see the implementation of the `Interpolate` routine in the `LS_Many_Ellipses` class for more info.

## Generate Mesh

Next, generate the mesh exactly how we did before:
```matlab
% setup up handle to interpolation routine
Interp_Handle = @(pt) LS.Interpolate(pt);

MG = MG.Get_Cut_Info(Interp_Handle);
[TRI, VTX] = MG.run_mex(Interp_Handle);
```
and plot it:
```matlab
trimesh(TRI,VTX(:,1),VTX(:,2));
axis equal;
grid on;
```

## How To Define Other Surfaces

To define your own domain, simply copy (and rename) the directory 
```matlab
./FELICITY/Static_Codes/Isosurface_Meshing/LevelSets_2D/@LS_Many_Ellipses
```
to another name of your choosing.  Note that `LS_Many_Ellipses` is sub-classed from the class Abstract_LevelSet, which contains internal variables `Param` and `Grid` for storing the data that defines the bounding surface.  Then, you just need to rewrite the `Interpolate` routine to match your domain.  Look at the unit tests in
```matlab
./FELICITY/Static_Codes/Isosurface_Meshing/Unit_Test/
```
for more info and examples.

*Note* you can mesh domains described by a polygon by using MATLAB's `inpolygon` routine to indicate whether a point is inside or outside.  However, you can only use bisection to compute the cut points.  Also note that the output mesh from the TIGER algorithm will not have the same vertices as the polygonal curve used to describe it.

# Simple Example (3-D)

Mesh a sphere.

## Initialization

Initializing the mesh generator is basically the same:
```matlab
% create mesher object
Cube_Dim = [0, 1];
Num_BCC_Points = 25;
Use_Newton = true;
TOL = 1e-12; % tolerance to use in computing the cut points
% BCC mesh of the unit cube [0,1] x [0,1] x [0,1]
MG = Mesher3Dmex(Cube_Dim,Num_BCC_Points,Use_Newton,TOL);
```

Here we use Newton's method for approximating the cut points.  Note, this means you must provide the level set function value, as well as its gradient, in the `Interpolate` routine.

The tolerance is set to 1E-12 and refers to the following convergence criteria:
```matlab
|f(p)| < TOL = 1E-12.
```
where `f` is the level set function value.  Note: we always assume we are meshing the volume bounded by the *zero* level set.

## Define Iso-surface

We use a FELICITY provided MATLAB class to define the sphere:
```matlab
LS = LS_Sphere();
LS.Param.cx   = 0.5;
LS.Param.cy   = 0.5;
LS.Param.cz   = 0.5;
LS.Param.rad  = 0.1;
LS.Param.sign = 1;
```
The sphere has radius 0.1, and it is centered at (0.5,0.5,0.5).

## Generate Mesh

Generating the mesh is the same as before:
```matlab
% setup up handle to interpolation routine
Interp_Handle = @(pt) LS.Interpolate(pt);

MG = MG.Get_Cut_Info(Interp_Handle);
[TET, VTX] = MG.run_mex(Interp_Handle);
```

This outputs the tetrahedron connectivity data `TET` (Mx4 matrix, where M is the number of tetrahedrons), and the vertex coordinates `VTX` (Nx3 matrix, where N is the number of vertices).

Now you can plot the mesh with the MATLAB command tetramesh:
```matlab
tetramesh(TET,VTX);
axis equal;
grid on;
```
Note: you can also use the FELICITY class `MeshTetrahedron` to plot it and make other manipulations.  See the tutorial [Mesh Classes](../wiki/Tutorial_Meshes_1) for more info.

# Volume Data (3-D)

Mesh volumetric data.  This is probably the most common use for the TIGER algorithm.

## Initialization

Initialize the mesh generator:
```matlab
% create mesher object
Cube_Dim = [0, 1];
Num_BCC_Points = 25;
Use_Newton = false;
TOL = 1e-2; % tolerance to use in computing the cut points
% BCC mesh of the unit cube [0,1] x [0,1] x [0,1]
MG = Mesher3Dmex(Cube_Dim,Num_BCC_Points,Use_Newton,TOL);
```

We only use bisection here.

## Define Iso-surface

We use a FELICITY class to store the volumetric data:
```matlab
LS = LS_Vol();
% define a 3-D cartesian grid for sampling the volume data:
s_vec = linspace(-0.1,1.1,101)';
[X, Y, Z] = meshgrid(s_vec,s_vec,s_vec);
LS.Grid.s_vec = s_vec;

% define grid data
R = 0.1;
C0 = 0.4;
V_pert = sin(2*pi*X).*sin(2*pi*Y).*cos(4*pi*Z);
LS.Grid.V = R - sqrt((X - 0.5).^2 + (Y - 0.5).^2 + (Z - 0.5).^2) + C0 * V_pert;
```
In other words, we defined volumetric data whose zero level set is a sphere, but then we added a _perturbation_ to that.

## Generate Mesh

Generating the mesh:
```matlab
% setup up handle to interpolation routine
Interp_Handle = @(pt) LS.Interpolate(pt);

MG = MG.Get_Cut_Info(Interp_Handle);
[TET, VTX] = MG.run_mex(Interp_Handle);
```
and plot it:
```matlab
tetramesh(TET,VTX);
axis equal;
grid on;
```

# Conclusion

For 2-D and 3-D, the definition of the iso-surface (i.e. zero level set) is made in the Interpolate routine.  For an example, look in the directory:
```matlab
./FELICITY/Static_Codes/Isosurface_Meshing/LevelSets_3D/@LS_Vol
```

See the next tutorial [Mesh Generation With TIGER: Part 2](../wiki/Mesh_Generation_with_TIGER_2) on how to generate bulk meshes from polygons and surface triangulations.