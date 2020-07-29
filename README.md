# The book: *“Sequential Software for the Robust Multigrid Technique”* by S.&nbsp;I.&nbsp;Martynenko

[![Fortran CI](https://github.com/simartynenko/Robust_Multigrid_Technique_2020/workflows/Fortran%20CI/badge.svg)](https://github.com/simartynenko/Robust_Multigrid_Technique_2020/actions?query=workflow%3A%22Fortran+CI%22)

The repository contains examples and programs from the book 
*“Sequential Software for the Robust Multigrid Technique”* by S.&nbsp;I.&nbsp;Martynenko, Triumph, Moscow, 2020.
Complete text of the book in Russian is available via the [link](../blob/master/book.pdf). Abstract in English, theoretical background and description of the source code files are provided below.

## Abstract of the book

This book covers sequential software for the __Robust Multigrid Technique__ (RMT) and presents basic concepts of modern numerical methods for mathematical modeling of physical and chemical processes in the computational continuum mechanics (thermal conductivity, chemical hydrodynamics, convective heat transfer, electrodynamics, etc.). FORTRAN language subroutines for solving the boundary value problems of the computational continuum mechanics by point and block Seidel method (Vanka-type smoother) are given in Chapters 1-4. The multigrid module `RMT_3D_2020` is described in detail in Chapter 2 of this book. The module contains subroutines for implementation of the problem-independent components of the Robust Multigrid Technique for solving the three-dimensional boundary value problems. Examples of the multigrid module `RMT_3D_2020` application in static and dynamic cycles are given in Chapters 3 and 4. Chapter “Short history of RMT” represents the historical development of the Robust Multigrid Technique.

Potential readers are graduate students and researchers working in applied and numerical mathematics as well as multigrid practitioners and software programmers for modeling physical and chemical processes in aviation and space industries, power engineering, chemical technology and other branches of mechanical engineering.

## Mathematical background of the Robust Multigrid Technique

All the required mathematical background is provided in the book:
S.&nbsp;I.&nbsp;Martynenko *“The Robust Multigrid Technique: For Black-Box Software”*, de Gruyter, Berlin, 2017 (https://www.degruyter.com/view/title/527481).

This book presents a detailed description of a robust pseudomultigrid algorithm for solving (initial-)boundary value problems on structured grids in a black-box manner. To overcome the problem of robustness, the presented __Robust Multigrid Technique__ (RMT) is based on the application of the essential multigrid principle in a single grid algorithm. It results in an extremely simple, very robust and highly parallel solver with close-to-optimal algorithmic complexity and the least number of problem-dependent components. Topics covered include an introduction to the mathematical principles of multigrid methods, a detailed description of RMT, results of convergence analysis and complexity, possible expansion on unstructured grids, numerical experiments and a brief description of multigrid software, parallel RMT and estimations of speed-up and efficiency of the parallel multigrid algorithms, and finally applications of RMT for the numerical solution of the incompressible Navier-Stokes equations.

## List of modules

* Multigrid modules `RMT_3D_2020` (§ 1, § 3 chapter 2) includes subprograms for main components RMT for solving 3D boundary value problems.
*The multigrid module is available as a static library only for Intel Fortran Compiler under Windows and GFortran under Linux operation systems.*

* Module `RMT_3D_UFG` (§ 2 chapter II, p. 70) includes the subprogram `U_Finest_Grid` for the uniform grid generation using (2.1)–(2.4).

* Module `RMT_3D_NFG` (§ 6 chapter II, p. 93) includes subprograms `Finest_Grid` and `CELLS` for the nonuniform grid generation. Finite volume faces are located between the vertices. 

* Module `RMT_3D_NFG2` (§ 6 chapter II, p. 97) includes subprograms `Finest_Grid` and `CELLS` for the nonuniform grid generation. Vertices are located between the finite volume faces.

* Module `RMT_3D_TIMES` includes a basic time measurement function. 

## List of example programs

1. Program `RMT_3D_2020_Example_01.f90` (§ 2 chapter I, p. 45): solution of discrete Poisson equation (1.6) with the boundary conditions (1.7) by point Seidel method. 

1. Program `RMT_3D_2020_Example_02_M.f90` (§ 3 chapter I, p. 54): solution of discrete Poisson 
equation (1.6) with the boundary conditions (1.7) by block Seidel method. 
Used modules: `RMT_3D_2020_Example_02_1`, `RMT_3D_2020_Example_02_2` (p. 60), `RMT_3D_TIMES`.

1. Program `RMT_3D_2020_Example_03.f90` (§ 2 chapter II, p. 75): formation of a multigrid structure generated by a uniform grid without information files.
 
1. Program `RMT_3D_2020_Example_04.f90` (§ 4 chapter II, p. 79): formation of a multigrid structure generated by a uniform grid with information files.

1. Program `RMT_3D_2020_Example_05.f90` (§ 5 chapter II, p. 82): example of calling the subroutine Index_Mapping for computing the index mapping. 

1. Program `RMT_3D_2020_Example_06.f90` (§ 5 chapter II, p. 88): approximation of derivatives on a multigrid structure.
Subroutine CELLS (§ 6 chapter II, p. 92):  generation a non-uniform grids with the same grid aspect ratio. 

1. Program `RMT_3D_2020_Example_07.f90` (§ 6 chapter II, стр. 95): formation of a multigrid structure generated by a non-uniform grids with the same grid aspect ratio. 
Used modules: `RMT_3D_NFG`.

1. Program `RMT_3D_2020_Example_08.f90` (§ 6 chapter II, p. 99): formation of a multigrid structure generated by a non-uniform grids with the same grid aspect ratio. 
Used modules: `RMT_3D_NFG2`.

1. Program `RMT_3D_2020_Example_09_M.f90` (§ 2 chapter III, p. 107): approximation of integral (3.4) on the multigrid structure generated by a uniform grid.
Used modules: `RMT_3D_2020_Example_09_1.f90` and `RMT_3D_2020_Example_09_2.f90`.

1. Program `RMT_3D_2020_Example_10_M.f90` (§ 2 chapter III, p. 111): approximation of integral (3.4) on the multigrid structure generated by a non-uniform grid.
Used modules: `RMT_3D_2020_Example_09_1`, `RMT_3D_2020_Example_09_2`, `RMT_3D_2020`, `RMT_3D_NFG`.

1. Program `RMT_3D_2020_Example_11_M.f90` (§ 2 chapter III, p. 114): approximation of integral (3.4) on the multigrid structure generated by a non-uniform grid.
Used modules: `RMT_3D_2020_Example_09_1`, `RMT_3D_2020_Example_09_2`, `RMT_3D_2020`, `RMT_3D_NFG2`.

1. Program `RMT_3D_2020_Example_12_M.f90` (§ 4 chapter III, p. 122): numerical solution of the Dirichlet problem (3.9) on a multigrid structure generated by a uniform grid. Smoother is point Seidel method.
Used modules: `RMT_3D_2020_Example_12_1`, `RMT_3D_2020_Example_12_2`, `RMT_3D_2020`, `RMT_3D_NFG2`, `RMT_3D_TIMES`.

1. Program `RMT_3D_2020_Example_13_M.f90` (§ 4 chapter III, p. 140): numerical solution of the Dirichlet problem (3.9) on a multigrid structure generated by a uniform grid. Smoother is block Seidel method.
Used modules: `RMT_3D_2020_Example_13_1`, `RMT_3D_2020_Example_13_2`, `RMT_3D_2020`, `RMT_3D_UFG`, `RMT_3D_TIMES`.

1. Program `RMT_3D_2020_Example_14.f90` (§ 4 chapter II, p. 150): index mapping in the dynamic cycle.
Used modules: `RMT_3D_UFG`.

1. Program `RMT_3D_2020_Example_15_M.f90` (§ 3 chapter IV, p. 153): approximation of integrals on the multigrid structure generated by the dynamic grids of the first levels.
Used modules: `RMT_3D_2020_Example_15_1`, `RMT_3D_2020_Example_15_2`, `RMT_3D_2020`, `RMT_3D_UFG`.

1. Program `RMT_3D_2020_Example_16_M.f90` (§ 4 chapter IV, p. 153): first multigrid iteration of the dynamic cycle.
Used modules: `RMT_3D_2020_Example_16_1`, `RMT_3D_2020_Example_16_2`, `RMT_3D_2020`, `RMT_3D_UFG`, `RMT_3D_TIMES`.
