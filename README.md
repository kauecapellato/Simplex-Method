\begin{document}
# The Simplex Method

This repository contains 2 folders: The first one is the Phase II Simplex and the second one is the Phases I and II Simplex. For each one, we'll briefly describe the files and how each function works.

The Simplex algorithm implementation was inspired by [[1]]().

--------------
## Summary


1. [Phase II Simplex](#1)
    1. Naive Implementation(#1.1)
    2. Revised Implementation(#1.2)
    3. Full Tableau Implementation(#1.3)
2. [Phases I and II Simplex](#2)

--------------

<a name="1"></a>

## 1\. Phase II Simplex - 3 Implementations:
Here we implemented the Simplex Method using its three most known versions: Naive, Revised and Tableau.

Each function returns two values: **ind**, that represents whether the problem has optimal solution (**ind = 0**) or the problem is unbounded (**ind = -1**); and **v**, that represents the optimal solution if it exists, or the direction vector in which the objective funcion goes to ![](https://latex.codecogs.com/gif.latex?-\infty).

<a name="1.1"></a>
### 1.1\. Naive Implementation
For the naive version, the input parameters are: 

```
[ind v] = simplex_ing(A,b,c,m,n,x,indB)
```

In which:
\\

--------------

<a name="2"></a>
### 2\. Phases I and II Simplex

--------------

# References
[1] Bertsimas D., Tsitsiklis N. J. Introduction to Linear Optmization (1997). Athena Scientific, Belmont, Massachusetts. {#ref1}
