/** \mainpage Shape Analysis for Union of Balls 

\section Introduction

Consider a set of balls, where balls are allowed to intersect. Analyzing spatial 
properties of the union of balls in this set is not an easy task, because of the
intersections. This library offers computing the volume of the union. 

The library employs weighted alpha shapes from CGAL and is based on 
inclusion-exclusion principle in order to solve this task effectively. The code 
is based on the documentation paper "Measuring Space Filling Diagrams and Voids" 
written by Herbert Edelsbrunner and Ping Fu in 1994.

*/


/** @example ExampleUsage.cpp 
This is a simple example of how to compute the volume of a set of balls.
*/