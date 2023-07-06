#Matrix.h
## Introduction

In this project you will implement your own library for processing numerical matrices in the C programming language. Matrices are one of the basic data structures in programming, e.g. they are used to represent table values, for computational tasks and neural networks. As part of the project you will learn more about matrices and solidify knowledge of structured programming.
--
##Task
Implement basic operations with matrices (partially described (#matrix-operations)): create_matrix (creation), remove_matrix (cleaning and destruction), eq_matrix (comparison), sum_matrix (addition), sub_matrix (subtraction), mult_matrix (multiplication), mult_number (multiplication by number), transpose (transpose), determinant (calculation of determinant), calc_complements (calculation of matrix of algebraic complements), inverse_matrix (finding inverse of the matrix).
--
- The library must be developed in C language of C11 standard using gcc compiler
- Do not use outdated and legacy language constructions and library functions. Pay attention to the legacy and obsolete marks in the official documentation on the language and the libraries used. Use the POSIX.1-2017 standard.
- Make it as a static library (with the s21_matrix.h header file)
- The library must be developed according to the principles of structured programming;
- Prepare full coverage of library functions code with unit-tests using the Check library
- Unit tests must cover at least 80% of each function (checked using gcov)  
- Provide a Makefile for building the library and tests (with targets all, clean, test, s21_matrix.a, gcov_report)
- The gcov_report target should generate a gcov report in the form of an html page. Unit tests must be run with gcov flags to do this 
- The matrix must be implemented as the structure described [above](#matrix-structure-in-c-language)
- Verifiable accuracy of the fractional part is up to 6 decimal places
