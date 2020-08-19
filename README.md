# OR2

## Operational Research Course Project
This is the project of the Operational Research 2 course, hold by professor M. Fischetti in the Department of Engineering (DEI), University of Padua. The purpose of the project is to introduce the theoretical concept regarding the Traveling Salesman Problem (TSP) as done during the course, to describe our group implemented methods and present some comparison between each other in term of execution time and solution cost.

The theoretical part of the course is divided in two big block: in the first part are introduced the exact methods to resolve the TSP, in particular some optimization of the Simplex Algorithm and its implementation with a licensed software: CPLEX; in the second part are presented different heuristics which goal is to resolve larger instances of the TSP in acceptable time at the expenses of solution cost.

The most interesting part is the application of Mixed Integer Programming methods (and the CPLEX software) to resolve a specific Linear Programming problem (TSP) and the comparison of exact methods results with heuristics in terms of execution time and solution cost.

## Compile

### Dependency

The program use IBM ILOG CPLEX, a licensed program, which can be downloaded following the [IBM instructions](https://www.ibm.com/products/ilog-cplex-optimization-studio). The program has been tested only for version 12.10.

### Linux

```
make mainTSP
```
In makefile can be checked and modify the parameter to compile and execute `mainTSP`.

### Windows

## Run

The parameter are passed like pair of key-value in the format `-key value`.
H
