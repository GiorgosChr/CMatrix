<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [QR Method for Finding Eigenvectors and Eigenvalues](#qr-method-for-finding-eigenvectors-and-eigenvalues)
  - [1. Introduction](#1-introduction)
  - [2. QR Decomposition](#2-qr-decomposition)
  - [3. QR Algorithm](#3-qr-algorithm)
  - [4. Finding Eigenvalues](#4-finding-eigenvalues)
  - [5. Finding Eigenvectors](#5-finding-eigenvectors)
  - [6. Example](#6-example)
    - [Iteration 1:](#iteration-1)
    - [Iteration 2:$$](#iteration-2)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# QR Method for Finding Eigenvectors and Eigenvalues

The QR method is an iterative numerical algorithm used to find the eigenvalues and eigenvectors of a square matrix. It is based on the concept of QR decomposition and can be used to solve the eigenvalue problem for both symmetric and nonsymmetric matrices.

## 1. Introduction

The eigenvalue problem for a matrix A can be written as:

$$A \cdot x = \lambda \cdot x$$

where $A$ is the matrix, $x$ is the eigenvector, and $\lambda$ is the corresponding eigenvalue. The QR method provides an iterative approach to finding these eigenvalues and eigenvectors.

## 2. QR Decomposition

The first step in the QR method is to decompose the matrix $A$ into the product of two matrices: $Q$ and $$R$. The QR decomposition can be expressed as:

$$A = Q \cdot R$$

where $Q$ is an orthogonal matrix and $R$ is an upper triangular matrix.

## 3. QR Algorithm

The QR algorithm is an iterative process that repeatedly applies the QR decomposition to the matrix $A$. The steps of the algorithm are as follows:

1. Start with an initial matrix $A_0$, which is the original matrix $A$.
2. Compute the QR decomposition: $A_k = Q_k \cdot R_k$.
3. Update the matrix $A$: $A_{k+1} = R_k \cdot Q_k$.
4. Repeat steps 2 and 3 until convergence is achieved.

The convergence of the QR algorithm occurs when the off-diagonal elements of $A_k$ become sufficiently small.

## 4. Finding Eigenvalues

To find the eigenvalues using the QR method, we iterate the QR algorithm until convergence. The eigenvalues can be obtained from the diagonal elements of the final matrix $A_k$.

## 5. Finding Eigenvectors

Once the eigenvalues are found, we can compute the eigenvectors corresponding to each eigenvalue. The eigenvectors can be obtained by solving the system of equations:

$$(A - \lambda \cdot I) \cdot x = 0$$

where $A$ is the matrix, $\lambda$ is the eigenvalue, $I$ is the identity matrix, and $x$ is the eigenvector.

## 6. Example

Let's consider a square matrix $A$:

$$A = \begin{bmatrix} 2 & -1 \\ 4 & 3 \end{bmatrix}$$

To find the eigenvalues and eigenvectors using the QR method, we follow these steps:

1. Initialize $A_0 = A$.
2. Compute the QR decomposition: $A_0 = Q_0 \cdot R_0$.
3. Update $A$: $A_1 = R_0 \cdot Q_0$.
4. Repeat steps 2 and 3 until convergence.

Let's perform the iterations:

### Iteration 1:

$$A_0 = \begin{bmatrix} 2 & -1 \\ 4 & 3 \end{bmatrix}$$

QR decomposition:

$$Q_0 = \begin{bmatrix} 0.447 & -0.894 \\ 0.894 & 0.447 \end{bmatrix}$$
$$R_0 = \begin{bmatrix} 4.472 & 1.341 \\ 0 & 2.236 \end{bmatrix}$$

Updated matrix:

$$A_1 = R_0 \cdot Q_0 = \begin{bmatrix} 2.828 & 2.828 \\ 0 & -0.447 \end{bmatrix}$$

### Iteration 2:$$

$$A_1 = \begin{bmatrix} 2.828 & 2.828 \\ 0 & -0.447 \end{bmatrix}$$

QR decomposition:

$$Q_1 = \begin{bmatrix} 0.997 & -0.077 \\ 0.077 & 0.997 \end{bmatrix}$$
$$R_1 = \begin{bmatrix} 2.828 & 2.828 \\ 0 & -0.447 \end{bmatrix}$$

Updated matrix:

$$A_2 = R_1 \cdot Q_1 = \begin{bmatrix} 2.828 & 0 \\ 0 & -0.447 \end{bmatrix}$$

After convergence, the final matrix $A_k$ is:

$$A_k = \begin{bmatrix} 2.828 & 0 \\ 0 & -0.447 \end{bmatrix}$$

The eigenvalues are the diagonal elements of $A_k$:

$$\lambda_1 = 2.828,$$
$$\lambda_2 = -0.447.$$

To find the eigenvectors, we solve the systems of equations:

For $\lambda_1 = 2.828$:

$$(A - \lambda_1 \cdot I) \cdot x = 0$$

$$\begin{bmatrix} 2 - 2.828 & -1 \\ 4 & 3 - 2.828 \end{bmatrix} \cdot x = 0$$

$$\begin{bmatrix} -0.828 & -1 \\ 4 & 0.172 \end{bmatrix} \cdot x = 0$$

Solving the system, we find the eigenvector $$x_1 = \begin{bmatrix} 0.702 \\ 0.712 \end{bmatrix}$$

For $$\lambda_2 = -0.447$$

$$(A - \lambda_2 \cdot I) \cdot x = 0$$

$$\begin{bmatrix} 2 - (-0.447) & -1 \\ 4 & 3 - (-0.447) \end{bmatrix} \cdot x = 0$$

$$\begin{bmatrix} 2.447 & -1 \\ 4 & 3.447 \end{bmatrix} \cdot x = 0$$

Solving the system, we find the eigenvector $$x_2 = \begin{bmatrix} 0.684 \\ -0.730 \end{bmatrix}$$

Therefore, the eigenvalues and eigenvectors of matrix $A$ are:

$$\lambda_1 = 2.828, \;x_1 = \begin{bmatrix} 0.702 \\ 0.712 \end{bmatrix}$$

$$\lambda_2 = -0.447,\; x_2 = \begin{bmatrix} 0.684 \\ -0.730 \end{bmatrix}$$
