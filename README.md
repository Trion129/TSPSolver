# TSP Solver
This is a TSP Solver I created for solving the assignment on Discrete optimizations course.
It has an attempted approach with Simulated annealing with 2-swap.

The other approach uses Lin Kernighan Algorithm.
It also uses a splay tree based data structure to speed up the search to bring it as close as possible to a O(NlogN) solution.

It is able to get a satisfactory result on large number of points in euclidean space.
e.g. solves problem with 80k cities within 3 seconds.
I was mostly focusing on the solving speed and building a scalable solver that might have some inaccuracy.

### References:

#### Splaytree data structure:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.49.570&rep=rep1&type=pdf

#### Lin-Kernighan heuristic:
http://akira.ruc.dk/~keld/research/LKH/LKH-2.0/DOC/LKH_REPORT.pdf