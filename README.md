# The Simplex Method

This repository contains 2 folders: The first one is the Phase II Simplex and the second one is the Phases I and II Simplex. For each one, we'll briefly describe the files and how each function works.

The Simplex algorithm implementation was inspired by [[1]]() and written in [Octave](https://www.gnu.org/software/octave/index) programming language.

--------------
## Summary


1. [Simplex: Phase II](#1)
    1. [Naive Implementation](#1.1)
    2. [Revised Implementation](#1.2)
    3. [Full Tableau Implementation](#1.3)
2. [Simplex: Phases I and II](#2)
    1. [Revised Implementation](#2.1)
    2. [Full Tableau Implementation](#2.2)

--------------

<a name="1"></a>

## 1\. Simplex: Phase II
Here we implemented the Phase II Simplex Method using its three most known versions: Naive, Revised and Tableau. The Phase II Simplex starts with a initial Basic Feasible Solution (BFS) **x**.

This folder contains 3 files: ```simplex_ing```, ```simplex_res``` and ```simplex_tab```.

Each function returns two values: **ind**, that represents whether the problem has optimal solution (**ind = 0**) or the problem is unbounded (**ind = -1**); and **v**, that represents the optimal solution if it exists, or the direction vector in which the objective funcion goes to -infinity.

<a name="1.1"></a>
### 1.1\. Naive Implementation
For the naive version, the main function is defined as follows: 

```
[ind v] = simplex_ing(A,b,c,m,n,x,indB)
```

![img1](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_1.PNG)

The problems treated for the 3 versions for the Phase II Simplex were obtained from [[2]]():

![img2](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_problem.PNG)

The input parameters are:

```
A = [4/5 6/5 1 0 0 0; 2 3 0 1 0 0; 2/3 2 0 0 1 0; 16/3 4 0 0 0 1];
c = [-40; -50; 0; 0; 0; 0];
b = [16; 30; 16; 64];
x = [0; 0; 16; 30; 16; 64];
m = size(A,1);
n = size(A,2);
indB = [3; 4; 5; 6];
```


--------------

<a name="1.2"></a>

### 1.2\. Revised Implementation
For the revised version, the main function is defined as follows: 

```
[ind v] = simplex_res(A,b,c,m,n,x,indB,Binv)
```
Where **Binv** stands for "B inverse", which is defined as:

![img3](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_2.PNG)


And for the following problem we have:

![img2](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_problem.PNG)

The input parameters are:

```
A = [4/5 6/5 1 0 0 0; 2 3 0 1 0 0; 2/3 2 0 0 1 0; 16/3 4 0 0 0 1];
c = [-40; -50; 0; 0; 0; 0];
b = [16; 30; 16; 64];
x = [0; 0; 16; 30; 16; 64];
m = size(A,1);
n = size(A,2);
indB = [3; 4; 5; 6];
Binv = eye(4); 
```


--------------

<a name="1.3"></a>

### 1.3\. Full Tableau Implementation

For the full tableau version, the main function is defined as follows: 

```
[ind v] = simplex_tab(indB,tableau)
```
Where the parameter **tableau** is as follows:

![img3](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_3.PNG)


And for the following problem we have:

![img2](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_problem.PNG)

The input parameters are:

```
indB = [3; 4; 5; 6];
tableau = [0 -40 -50 0 0 0 0;
           16 4/5 6/5 1 0 0 0;
           30 2 3 0 1 0 0;
           16 2/3 2 0 0 1 0;
           64 16/3 4 0 0 0 1];
```


--------------
<a name="2"></a>
## 2\. Phases I and II Simplex

In this implementation, the Phase I Simplex tries to find a BFS that is used as the parameter **x** in the Phase II Simplex. The last two versions are used: Revised and Full Tableau.

As an anticiclying rule, this implementation uses the Bland's rule, that selects the smallest subscript as pivot.

This folder contains 4 files: ```simplex_res``` and ```simplex_revisado``` for the two Phase Simplex revised version and ```simplex_tab``` and ```simplex_tableau``` for the two Phase Simplex tableau version.

<a name="2.1"></a>
### 2.1\. Revised Implementation

The function is defined as:

```
[ind x v] = simplex_res(A,b,c,m,n)
```

Where x is the initial BFS found in the end of Phase I (for problems with optimal solution or unbounded problems) and v is the same as described on Section [1](#1) and is used for both versions (revised and tableau). Where:

![img4](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_4.PNG)


For the following problem:

![img4](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_problem2.PNG)

We have: 

```
A = [1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 1];
b = [3; 2; 5; 1];
c = [1; 1; 1; 0];
m = size(A, 1);
n = size(A,2);
```

--------------

<a name="2.2"></a>
### 2.2\. Full Tableau Implementation

For the tableau version, there is an unique parameter which is:

![img5](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_5.PNG)

For the problem:

![img4](https://github.com/kauecapellato/TheSimplexMethod/blob/main/simplex_images/github_problem2.PNG)

The input is:

```
tableau = [0 1 1 1 0;
           3 1 2 3 0;
           2 -1 2 6 0;
           5 0 4 9 0;
           1 0 0 3 1];
```
--------------

# References
[1] Bertsimas D., Tsitsiklis N. J. Introduction to Linear Optmization (1997). Athena Scientific, Belmont, Massachusetts.

[2] Graves, S. The Simplex Method for Solving a Linear Program (2003). Massachusetts Institute of Technology. Pages 55 - 58.
