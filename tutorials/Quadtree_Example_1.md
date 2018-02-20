Quadtree Implementation
=======================

# Introduction

FELICITY provides implementations of binary trees, quadtrees, and octrees for storing points distributed in space.  This can be used to query nearest neighbors, and find closest points.

# Quadtree Example (2-D)

Use a quadtree to efficiently find k nearest neighbors.

## Initialization

After you have installed FELICITY, make sure you run `test_FELICITY.m` or `compile_static_codes.m`.  This will compile the C++ code that implements the search tree.

Let us store some random 2-D points in a quadtree.  Type the following at the MATLAB prompt (or put it in a script file):
```matlab
% create a set of random 2-D points inside the unit square [0,1] x [0,1]
NUM = 10000;
points = rand(NUM,2);
```

Next, create a Quadtree object to store these points.  This requires a bounding box that *contains* the points.
```matlab
% create Quadtree...
BB = [-0.001, 1.001, -0.001, 1.001];
Max_Tree_Levels = 32;
Bucket_Size = 20;
QT = mexQuadtree(points,BB);%,Max_Tree_Levels,Bucket_Size);
clear points;
```
Note: one can provide additional arguments specifying the max tree depth and the bucket size (max number of points to store in each leaf cell).

Note: `QT.Points` contains a copy of the matrix `points`.

This creates a MATLAB object that runs the quadtree and interfaces with a C++ class object.

## Query Nearest Neighbors

Next, for each point in a given set of query points, we want to find the k nearest neighbors amongst the points we stored in the quadtree.  Thus, define some random query points:
```matlab
QP_Num = 100;
query_points = rand(QP_Num,2);
```

Now compute the 2 nearest neighbors for each query point:
```matlab
NN = 2;
[QT_indices, QT_dist] = QT.kNN_Search(query_points,NN);
```
Here, the ith row of `QT_indices` and `QT_dist` corresponds to the ith row of `query_points`.  The number of columns of `QT_indices` and `QT_dist` is 2; the first column corresponds to the closest neighbor, the second column corresponds to the next closest neighbor.

`QT_indices` lists row indices (into `QT.Points`), i.e. this indicates which of the points in `QT.Points` are the nearest neighbors.

For example, `QT_indices(i,:)` gives the two nearest neighbors (in `QT.Points`) to the ith point in `query_points`.

`QT_dist` gives the euclidean distance to the nearest neighbors.

## Change The Point Positions

Suppose you want to change the coordinates of the points stored in the quadtree.  This could be because the points are moving in a simulation.  Here is how to do it.
```matlab
% generate new points
new_pts = rand(NUM,2);
% Update the tree with brand new points
QT = QT.Update_Tree(new_pts);
clear new_pts;
```
Note: `new_pts` must have the same `size` as `QT.Points`.

Now recompute the nearest neighbors:
```matlab
[QT_indices, QT_dist] = QT.kNN_Search(query_points,NN);
```

Make a plot of the 2 nearest neighbors to `query_points(1,:)`:
```matlab
figure;
plot(QT.Points(:,1),QT.Points(:,2),'k*');
hold on;
plot(query_points(1,1),query_points(1,2),'bp');
plot(QT.Points(QT_indices(1,:),1),QT.Points(QT_indices(1,:),2),'r.','MarkerSize',20);
hold off;
axis([0 1 0 1]);
axis equal;
```

## Final Cleanup

When you are done with the object, you *must* delete it:
```matlab
% delete the C++ object
delete(QT);
```

# Conclusion

See the unit tests in the directory:
```matlab
./FELICITY/Static_Codes/Search_Trees/Unit_Test/
```
for more information and examples.