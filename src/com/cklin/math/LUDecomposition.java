/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * This file is is modified by kkk on 2017.
 */

package com.cklin.math;

/**
 * Calculates the LUP-decomposition of a square matrix.
 * <p>The LUP-decomposition of a matrix A consists of three matrices L, U and
 * P that satisfy: P&times;A = L&times;U. L is lower triangular (with unit
 * diagonal terms), U is upper triangular and P is a permutation matrix. All
 * matrices are m&times;m.</p>
 * <p>As shown by the presence of the P matrix, this decomposition is
 * implemented using partial pivoting.</p>
 * <p>This class is based on the class with similar name from the
 * <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</p>
 * <ul>
 * <li>a {@link #getP() getP} method has been added,</li>
 * <li>the {@code det} method has been renamed as {@link #getDeterminant()
 * getDeterminant},</li>
 * <li>the {@code getDoublePivot} method has been removed (but the int based
 * {@link #getPivot() getPivot} method has been kept),</li>
 * <li>the {@code solve} and {@code isNonSingular} methods have been replaced
 * by a {@link #getSolver() getSolver} method and the equivalent methods
 * provided by the returned {@link DecompositionSolver}.</li>
 * </ul>
 *
 * @see <a href="http://mathworld.wolfram.com/LUDecomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/LU_decomposition">Wikipedia</a>
 * @since 2.0 (changed to concrete class in 3.0)
 */
public class LUDecomposition {
    /**
     * Default bound to determine effective singularity in LU decomposition.
     */
    private static final double DEFAULT_TOO_SMALL = 1e-16;
    /**
     * Default bound to determine effective singularity in LU decomposition.
     */
    private static final double WARING_TOO_SMALL = 1e-11;
    /**
     * Entries of LU decomposition.
     */
    private final double[][] lu;
    /**
     * Pivot permutation associated with LU decomposition.
     */
    private final int[] pivot;
    /**
     * Parity of the permutation associated with the LU decomposition.
     */
    private boolean even;
    /**
     * Singularity indicator.
     */
    private boolean singular;

    /**
     * Calculates the LU-decomposition of the given matrix.
     * This constructor uses 1e-11 as default value for the singularity
     * threshold.
     *
     * @param matrix Matrix to decompose.
     * @throws NonSquareMatrixException if matrix is not square.
     */
    public LUDecomposition(double[][] A) {
        this(A, DEFAULT_TOO_SMALL);
    }

    /**
     * Calculates the LU-decomposition of the given matrix.
     *
     * @param matrix               The matrix to decompose.
     * @param singularityThreshold threshold (based on partial row norm)
     *                             under which a matrix is considered singular
     * @throws NonSquareMatrixException if matrix is not square
     */
    public LUDecomposition(double[][] A, double singularityThreshold) {
        if (A.length != A[0].length)
            throw new IllegalArgumentException("A must be SquareMatrix");
        ;
        final int m = ArraysUtils.colCount(A);
        lu = ArraysUtils.copy(A);
        pivot = new int[m];

        // Initialize permutation array and parity
        for (int row = 0; row < m; row++) {
            pivot[row] = row;
        }
        even = true;
        singular = false;

        // Loop over columns
        for (int col = 0; col < m; col++) {

            // upper
            for (int row = 0; row < col; row++) {
                final double[] luRow = lu[row];
                double sum = luRow[col];
                for (int i = 0; i < row; i++) {
                    sum -= luRow[i] * lu[i][col];
                }
                luRow[col] = sum;
            }

            // lower
            int max = col; // permutation row
            double largest = Double.NEGATIVE_INFINITY;
            for (int row = col; row < m; row++) {
                final double[] luRow = lu[row];
                double sum = luRow[col];
                for (int i = 0; i < col; i++) {
                    sum -= luRow[i] * lu[i][col];
                }
                luRow[col] = sum;

                // maintain best permutation choice
                if (Math.abs(sum) > largest) {
                    largest = Math.abs(sum);
                    max = row;
                }
            }

            // Singularity check
            if (Math.abs(lu[max][col]) < singularityThreshold) {
                singular = true;
                return;
            }
            if (Math.abs(lu[max][col]) < WARING_TOO_SMALL) {
                System.out.println("Warning: Matrix is singular to working precision");
            }

            // Pivot if necessary
            if (max != col) {
                double tmp = 0;
                final double[] luMax = lu[max];
                final double[] luCol = lu[col];
                for (int i = 0; i < m; i++) {
                    tmp = luMax[i];
                    luMax[i] = luCol[i];
                    luCol[i] = tmp;
                }
                int temp = pivot[max];
                pivot[max] = pivot[col];
                pivot[col] = temp;
                even = !even;
            }

            // Divide the lower elements by the "winning" diagonal elt.
            final double luDiag = lu[col][col];
            for (int row = col + 1; row < m; row++) {
                lu[row][col] /= luDiag;
            }
        }
    }

    /**
     * Returns the matrix L of the decomposition.
     * <p>L is a lower-triangular matrix</p>
     *
     * @return the L matrix (or null if decomposed matrix is singular)
     */
    public double[][] getL() {
        final int m = pivot.length;
        double[][] L = ArraysUtils.zeros(m, m);
        for (int i = 0; i < m; ++i) {
            final double[] luI = lu[i];
            for (int j = 0; j < i; ++j) {
                L[i][j] = luI[j];
            }
            L[i][i] = 1.0;
        }
        return L;
    }

    /**
     * Returns the matrix U of the decomposition.
     * <p>U is an upper-triangular matrix</p>
     *
     * @return the U matrix (or null if decomposed matrix is singular)
     */
    public double[][] getU() {
        final int m = pivot.length;
        double[][] U = ArraysUtils.zeros(m, m);
        for (int i = 0; i < m; ++i) {
            final double[] luI = lu[i];
            for (int j = i; j < m; ++j) {
                U[i][j] = luI[j];
            }
        }
        return U;
    }

    /**
     * Returns the P rows permutation matrix.
     * <p>P is a sparse matrix with exactly one element set to 1.0 in
     * each row and each column, all other elements being set to 0.0.</p>
     * <p>The positions of the 1 elements are given by the {@link #getPivot()
     * pivot permutation vector}.</p>
     *
     * @return the P rows permutation matrix (or null if decomposed matrix is singular)
     * @see #getPivot()
     */
    public double[][] getP() {
        final int m = pivot.length;
        double[][] P = ArraysUtils.zeros(m, m);
        for (int i = 0; i < m; ++i)
            P[i][pivot[i]] = 1.0;
        return P;
    }

    /**
     * Returns the pivot permutation vector.
     *
     * @return the pivot permutation vector
     * @see #getP()
     */
    public int[] getPivot() {
        return pivot.clone();
    }

    /**
     * Return the determinant of the matrix
     *
     * @return determinant of the matrix
     */
    public double getDeterminant() {
        if (singular) {
            return 0;
        } else {
            final int m = pivot.length;
            double determinant = even ? 1 : -1;
            for (int i = 0; i < m; i++) {
                determinant *= lu[i][i];
            }
            return determinant;
        }
    }

    /**
     * Solve the linear equation A &times; X = B for matrices A.
     * <p>
     * The A matrix is implicit, it is provided by the underlying
     * decomposition algorithm.
     *
     * @param b right-hand side of the equation A &times; X = B
     * @return a vector X that minimizes the two norm of A &times; X - B
     * @throws org.apache.commons.math3.exception.DimensionMismatchException if the matrices dimensions do not match.
     * @throws SingularMatrixException                                       if the decomposed matrix is singular.
     */
    public double[] solve(double[] b) {
        final int m = pivot.length;
        if (b.length != m) {
            throw new IllegalArgumentException("Matrix dimensions must agree: " + m);
        }
        if (singular) {
            throw new RuntimeException("Matrix is singular.");
        }

        final double[] bp = new double[m];

        // Apply permutations to b
        for (int row = 0; row < m; row++) {
            bp[row] = b[pivot[row]];
        }

        // Solve LY = b
        for (int col = 0; col < m; col++) {
            final double bpCol = bp[col];
            for (int i = col + 1; i < m; i++) {
                bp[i] -= bpCol * lu[i][col];
            }
        }

        // Solve UX = Y
        for (int col = m - 1; col >= 0; col--) {
            bp[col] /= lu[col][col];
            final double bpCol = bp[col];
            for (int i = 0; i < col; i++) {
                bp[i] -= bpCol * lu[i][col];
            }
        }

        return bp;
    }

    /**
     * Solve the linear equation A &times; X = B for matrices A.
     * <p>
     * The A matrix is implicit, it is provided by the underlying
     * decomposition algorithm.
     *
     * @param b right-hand side of the equation A &times; X = B
     * @return a matrix X that minimizes the two norm of A &times; X - B
     * @throws org.apache.commons.math3.exception.DimensionMismatchException if the matrices dimensions do not match.
     * @throws SingularMatrixException                                       if the decomposed matrix is singular.
     */
    public double[][] solve(double[][] b) {

        final int m = pivot.length;
        if (ArraysUtils.rowCount(b) != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (singular) {
            throw new RuntimeException("Matrix is singular.");
        }

        final int nColB = ArraysUtils.colCount(b);

        // Apply permutations to b
        final double[][] bp = new double[m][nColB];
        for (int row = 0; row < m; row++) {
            final double[] bpRow = bp[row];
            final int pRow = pivot[row];
            for (int col = 0; col < nColB; col++) {
                bpRow[col] = b[pRow][col];
            }
        }

        // Solve LY = b
        for (int col = 0; col < m; col++) {
            final double[] bpCol = bp[col];
            for (int i = col + 1; i < m; i++) {
                final double[] bpI = bp[i];
                final double luICol = lu[i][col];
                for (int j = 0; j < nColB; j++) {
                    bpI[j] -= bpCol[j] * luICol;
                }
            }
        }

        // Solve UX = Y
        for (int col = m - 1; col >= 0; col--) {
            final double[] bpCol = bp[col];
            final double luDiag = lu[col][col];
            for (int j = 0; j < nColB; j++) {
                bpCol[j] /= luDiag;
            }
            for (int i = 0; i < col; i++) {
                final double[] bpI = bp[i];
                final double luICol = lu[i][col];
                for (int j = 0; j < nColB; j++) {
                    bpI[j] -= bpCol[j] * luICol;
                }
            }
        }
        return bp;
    }

    /**
     * Get the inverse of the decomposed matrix.
     *
     * @return the inverse matrix.
     * @throws SingularMatrixException if the decomposed matrix is singular.
     */
    public double[][] getInverse() {
        if (singular)
            throw new RuntimeException("Matrix is singular.");
        return solve(ArraysUtils.eye(pivot.length));
    }
}