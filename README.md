# Cayley-Hamilton Theorem Verifier

## Overview
This project implements an application that generates a random `3 × 3` matrix with entries from the finite field `Z_5`, computes its characteristic polynomial, determines its eigenvalues, and verifies the Cayley-Hamilton theorem. If the matrix has three distinct eigenvalues, the app also computes the corresponding eigenvectors.

## Features
- Generates a random `3 × 3` matrix with elements from `Z_5`
- Computes the characteristic polynomial of the matrix
- Finds the eigenvalues by evaluating the polynomial at `t = 0, 1, 2, 3, 4`
- Computes eigenvectors if there are three distinct eigenvalues
- Verifies the Cayley-Hamilton theorem

## Mathematical Background
The **Cayley-Hamilton theorem** states that every square matrix satisfies its own characteristic equation. Given a `3 × 3` matrix `A`, its characteristic polynomial is given by:
```
p(t) = det(A - tI)
```
where `I` is the identity matrix. The coefficients of the polynomial include:
- The determinant of `A`
- The trace of `A`
- The constant term `1` (up to sign)

The eigenvalues of `A` are the values of `t` in `Z_5` for which `p(t) = 0`. If the characteristic polynomial factors into distinct linear terms, we can compute eigenvectors accordingly.

## Implementation Details
- The program constructs the matrix with random entries from `Z_5`
- Computes `p(t)` and extracts its coefficients
- Finds eigenvalues by evaluating `p(t)` at `t = 0, 1, 2, 3, 4`
- If there are three distinct eigenvalues, it computes the corresponding eigenvectors
- Finally, it checks that `A` satisfies its characteristic equation
