Mesh classes in FELICITY
========================

FELICITY has classes that are convenient for _manipulating_ meshes. Once FELICITY is installed (and is in the MATLAB path!), type the following at the MATLAB prompt:

```matlab
>> Vtx = [0 0 0; 1 0 0; 0 1 1]; % 3-D coordinates
>> Tri = [1 2 3]; % triangle connectivity
>> Mesh = MeshTriangle(Tri,Vtx,'Omega');
```

This creates a `MeshTriangle` object. Now type the following and press "ENTER":

```matlab
>> Mesh
```

This causes MATLAB to display the members of the object. Also, try entering:

```matlab
>> methods(Mesh)
```

This will display the various methods available for that object. One important method is uniform refinement:

```matlab
>> Mesh = Mesh.Refine;
```

This refines the (initial) single triangle into four self-similar triangles. You can view it with the command:

```matlab
>> Mesh.Plot;
```

You can also refine the mesh *selectively* using Rivara's bisection algorithm. Define an array of _marked_ triangles and input that to the refinement method:

```matlab
>> Marked = [1; 3; 4];
>> Mesh = Mesh.Refine('bisection',Marked);
```

This will cause the first, third, and fourth triangles in the mesh to be bisected along their longest edge. The algorithm will then bisect any other triangles necessary to maintain a conforming mesh. 

For more information, see the *Basic Mesh Manipulation* chapter in the PDF manual.