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

import java.util.Arrays;


/**
 * Calculates the QR-decomposition of a matrix.
 * <p>The QR-decomposition of a matrix A consists of two matrices Q and R
 * that satisfy: A = QR, Q is orthogonal (Q<sup>T</sup>Q = I), and R is
 * upper triangular. If A is m&times;n, Q is m&times;m and R m&times;n.</p>
 * <p>This class compute the decomposition using Householder reflectors.</p>
 * <p>For efficiency purposes, the decomposition in packed form is transposed.
 * This allows inner loop to iterate inside rows, which is much more cache-efficient
 * in Java.</p>
 * <p>This class is based on the class with similar name from the
 * <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
 * following changes:</p>
 * <ul>
 * <li>a {@link #getQT() getQT} method has been added,</li>
 * <li>the {@code solve} and {@code  } methods have been replaced
 * by a {@link #getSolver() getSolver} method and the equivalent methods
 * provided by the returned {@link DecompositionSolver}.</li>
 * </ul>
 *
 * @see <a href="http://mathworld.wolfram.com/QRDecomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/QR_decomposition">Wikipedia</a>
 * @since 1.2 (changed to concrete class in 3.0)
 */
public class QRDecomposition {
    /**
     * A packed TRANSPOSED representation of the QR decomposition.
     * <p>The elements BELOW the diagonal are the elements of the UPPER triangular
     * matrix R, and the rows ABOVE the diagonal are the Householder reflector vectors
     * from which an explicit form of Q can be recomputed if desired.</p>
     */
    private double[][] qrt;
    /**
     * The diagonal elements of R.
     */
    private double[] rDiag;
    /**
     * Singularity threshold.
     */
    private final double threshold;

    /**
     * Calculates the QR-decomposition of the given matrix.
     * The singularity threshold defaults to zero.
     *
     * @param matrix The matrix to decompose.
     * @see #QRDecomposition(RealMatrix, double)
     */
    public QRDecomposition(double[][] matrix) {
        this(matrix, 0d);
    }

    /**
     * Calculates the QR-decomposition of the given matrix.
     *
     * @param matrix    The matrix to decompose.
     * @param threshold Singularity threshold.
     */
    public QRDecomposition(double[][] matrix,
                           double threshold) {
        this.threshold = threshold;

        final int m = matrix.length;
        final int n = matrix[0].length;
        qrt = ArraysUtils.transpose(matrix);
        rDiag = new double[Math.min(m, n)];

        decompose(qrt);
    }

    /**
     * Decompose matrix.
     *
     * @param matrix transposed matrix
     * @since 3.2
     */
    protected void decompose(double[][] matrix) {
        for (int minor = 0; minor < Math.min(matrix.length, matrix[0].length); minor++) {
            performHouseholderReflection(minor, matrix);
        }
    }

    /**
     * Perform Householder reflection for a minor A(minor, minor) of A.
     *
     * @param minor  minor index
     * @param matrix transposed matrix
     * @since 3.2
     */
    protected void performHouseholderReflection(int minor, double[][] matrix) {

        final double[] qrtMinor = matrix[minor];

		/*
         * Let x be the first column of the minor, and a^2 = |x|^2.
		 * x will be in the positions qr[minor][minor] through qr[m][minor].
		 * The first column of the transformed minor will be (a,0,0,..)'
		 * The sign of a is chosen to be opposite to the sign of the first
		 * component of x. Let's find a:
		 */
        double xNormSqr = 0;
        for (int row = minor; row < qrtMinor.length; row++) {
            final double c = qrtMinor[row];
            xNormSqr += c * c;
        }
        final double a = (qrtMinor[minor] > 0) ? -Math.sqrt(xNormSqr) : Math.sqrt(xNormSqr);
        rDiag[minor] = a;

        if (a != 0.0) {

			/*
			 * Calculate the normalized reflection vector v and transform
			 * the first column. We know the norm of v beforehand: v = x-ae
			 * so |v|^2 = <x-ae,x-ae> = <x,x>-2a<x,e>+a^2<e,e> =
			 * a^2+a^2-2a<x,e> = 2a*(a - <x,e>).
			 * Here <x, e> is now qr[minor][minor].
			 * v = x-ae is stored in the column at qr:
			 */
            qrtMinor[minor] -= a; // now |v|^2 = -2a*(qr[minor][minor])

			/*
			 * Transform the rest of the columns of the minor:
			 * They will be transformed by the matrix H = I-2vv'/|v|^2.
			 * If x is a column vector of the minor, then
			 * Hx = (I-2vv'/|v|^2)x = x-2vv'x/|v|^2 = x - 2<x,v>/|v|^2 v.
			 * Therefore the transformation is easily calculated by
			 * subtracting the column vector (2<x,v>/|v|^2)v from x.
			 *
			 * Let 2<x,v>/|v|^2 = alpha. From above we have
			 * |v|^2 = -2a*(qr[minor][minor]), so
			 * alpha = -<x,v>/(a*qr[minor][minor])
			 */
            for (int col = minor + 1; col < matrix.length; col++) {
                final double[] qrtCol = matrix[col];
                double alpha = 0;
                for (int row = minor; row < qrtCol.length; row++) {
                    alpha -= qrtCol[row] * qrtMinor[row];
                }
                alpha /= a * qrtMinor[minor];

                // Subtract the column vector alpha*v from x.
                for (int row = minor; row < qrtCol.length; row++) {
                    qrtCol[row] -= alpha * qrtMinor[row];
                }
            }
        }
    }


    /**
     * Returns the matrix R of the decomposition.
     * <p>R is an upper-triangular matrix</p>
     *
     * @return the R matrix
     */
    public double[][] getR() {
        // R is supposed to be m x n
        final int n = qrt.length;
        final int m = qrt[0].length;
        double[][] ra = new double[m][n];
        // copy the diagonal from rDiag and the upper triangle of qr
        for (int row = Math.min(m, n) - 1; row >= 0; row--) {
            ra[row][row] = rDiag[row];
            for (int col = row + 1; col < n; col++) {
                ra[row][col] = qrt[col][row];
            }
        }
        return ra;
    }

    /**
     * Returns the matrix Q of the decomposition.
     * <p>Q is an orthogonal matrix</p>
     *
     * @return the Q matrix
     */
    public double[][] getQ() {
        return ArraysUtils.transpose(getQT());
    }

    /**
     * Returns the transpose of the matrix Q of the decomposition.
     * <p>Q is an orthogonal matrix</p>
     *
     * @return the transpose of the Q matrix, Q<sup>T</sup>
     */
    public double[][] getQT() {
        // QT is supposed to be m x m
        final int n = qrt.length;
        final int m = qrt[0].length;
        double[][] qta = new double[m][m];

		/*
		 * Q = Q1 Q2 ... Q_m, so Q is formed by first constructing Q_m and then
		 * applying the Householder transformations Q_(m-1),Q_(m-2),...,Q1 in
		 * succession to the result
		 */
        for (int minor = m - 1; minor >= Math.min(m, n); minor--) {
            qta[minor][minor] = 1.0d;
        }

        for (int minor = Math.min(m, n) - 1; minor >= 0; minor--) {
            final double[] qrtMinor = qrt[minor];
            qta[minor][minor] = 1.0d;
            if (qrtMinor[minor] != 0.0) {
                for (int col = minor; col < m; col++) {
                    double alpha = 0;
                    for (int row = minor; row < m; row++) {
                        alpha -= qta[col][row] * qrtMinor[row];
                    }
                    alpha /= rDiag[minor] * qrtMinor[minor];

                    for (int row = minor; row < m; row++) {
                        qta[col][row] += -alpha * qrtMinor[row];
                    }
                }
            }
        }
        return qta;
    }

    /**
     * Returns the Householder reflector vectors.
     * <p>H is a lower trapezoidal matrix whose columns represent
     * each successive Householder reflector vector. This matrix is used
     * to compute Q.</p>
     *
     * @return a matrix containing the Householder reflector vectors
     */
    public double[][] getH() {
        final int n = qrt.length;
        final int m = qrt[0].length;
        double[][] ha = new double[m][n];
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < Math.min(i + 1, n); ++j) {
                ha[i][j] = qrt[j][i] / -rDiag[j];
            }
        }
        return ha;
    }

    /**
     * Check if the decomposed matrix is non-singular.
     *
     * @return true if the decomposed matrix is non-singular.
     */
    public boolean isNonSingular() {
        for (double diag : rDiag) {
            if (Math.abs(diag) <= threshold) {
                return false;
            }
        }
        return true;
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
        final int n = qrt.length;
        final int m = qrt[0].length;
        if (b.length != m) {
            throw new IllegalArgumentException("Matrix dimensions must agree: " + m);
        }
        if (!isNonSingular()) {
            throw new RuntimeException("Matrix is singular.");
        }

        final double[] x = new double[n];
        final double[] y = b.clone();

        // apply Householder transforms to solve Q.y = b
        for (int minor = 0; minor < Math.min(m, n); minor++) {

            final double[] qrtMinor = qrt[minor];
            double dotProduct = 0;
            for (int row = minor; row < m; row++) {
                dotProduct += y[row] * qrtMinor[row];
            }
            dotProduct /= rDiag[minor] * qrtMinor[minor];

            for (int row = minor; row < m; row++) {
                y[row] += dotProduct * qrtMinor[row];
            }
        }

        // solve triangular system R.x = y
        for (int row = rDiag.length - 1; row >= 0; --row) {
            y[row] /= rDiag[row];
            final double yRow = y[row];
            final double[] qrtRow = qrt[row];
            x[row] = yRow;
            for (int i = 0; i < row; i++) {
                y[i] -= yRow * qrtRow[i];
            }
        }

        return x;
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
        final int n = qrt.length;
        final int m = qrt[0].length;
        if (b.length != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!isNonSingular()) {
            throw new RuntimeException("Matrix is singular.");
        }

        final int columns = b[0].length;
        final int blockSize = BLOCK_SIZE;
        final int cBlocks = (columns + blockSize - 1) / blockSize;
        final double[][] xBlocks = createBlocksLayout(n, columns);
        final double[] alpha = new double[blockSize];

        for (int kBlock = 0; kBlock < cBlocks; ++kBlock) {
            final int kStart = kBlock * blockSize;
            final int kEnd = Math.min(kStart + blockSize, columns);
            final int kWidth = kEnd - kStart;

            // get the right hand side vector
            double[][] y = ArraysUtils.getCols(b, ArraysUtils.arithmetic(kStart, 1, kEnd - 1));
            //b.copySubMatrix(0, m - 1, kStart, kEnd - 1, y);

            // apply Householder transforms to solve Q.y = b
            for (int minor = 0; minor < Math.min(m, n); minor++) {
                final double[] qrtMinor = qrt[minor];
                final double factor = 1.0 / (rDiag[minor] * qrtMinor[minor]);

                Arrays.fill(alpha, 0, kWidth, 0.0);
                for (int row = minor; row < m; ++row) {
                    final double d = qrtMinor[row];
                    final double[] yRow = y[row];
                    for (int k = 0; k < kWidth; ++k) {
                        alpha[k] += d * yRow[k];
                    }
                }
                for (int k = 0; k < kWidth; ++k) {
                    alpha[k] *= factor;
                }

                for (int row = minor; row < m; ++row) {
                    final double d = qrtMinor[row];
                    final double[] yRow = y[row];
                    for (int k = 0; k < kWidth; ++k) {
                        yRow[k] += alpha[k] * d;
                    }
                }
            }

            // solve triangular system R.x = y
            for (int j = rDiag.length - 1; j >= 0; --j) {
                final int jBlock = j / blockSize;
                final int jStart = jBlock * blockSize;
                final double factor = 1.0 / rDiag[j];
                final double[] yJ = y[j];
                final double[] xBlock = xBlocks[jBlock * cBlocks + kBlock];
                int index = (j - jStart) * kWidth;
                for (int k = 0; k < kWidth; ++k) {
                    yJ[k] *= factor;
                    xBlock[index++] = yJ[k];
                }

                final double[] qrtJ = qrt[j];
                for (int i = 0; i < j; ++i) {
                    final double rIJ = qrtJ[i];
                    final double[] yI = y[i];
                    for (int k = 0; k < kWidth; ++k) {
                        yI[k] -= yJ[k] * rIJ;
                    }
                }
            }
        }
        return BlocksToArrays(n, columns, xBlocks);
    }

    //Commons Math's BlockMatrix
    private final int BLOCK_SIZE = 52;

    private double[][] BlocksToArrays(int rows, int columns, double[][] blocks) {
        // number of blocks
        final int blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
        final int blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;
        final double[][] data = new double[rows][columns];
        final int lastColumns = columns - (blockColumns - 1) * BLOCK_SIZE;

        for (int iBlock = 0; iBlock < blockRows; ++iBlock) {
            final int pStart = iBlock * BLOCK_SIZE;
            final int pEnd = Math.min(pStart + BLOCK_SIZE, rows);
            int regularPos = 0;
            int lastPos = 0;
            for (int p = pStart; p < pEnd; ++p) {
                final double[] dataP = data[p];
                int blockIndex = iBlock * blockColumns;
                int dataPos = 0;
                for (int jBlock = 0; jBlock < blockColumns - 1; ++jBlock) {
                    System.arraycopy(blocks[blockIndex++], regularPos, dataP, dataPos, BLOCK_SIZE);
                    dataPos += BLOCK_SIZE;
                }
                System.arraycopy(blocks[blockIndex], lastPos, dataP, dataPos, lastColumns);
                regularPos += BLOCK_SIZE;
                lastPos += lastColumns;
            }
        }

        return data;
    }

    private double[][] createBlocksLayout(final int rows, final int columns) {
        final int blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
        final int blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE;

        final double[][] blocks = new double[blockRows * blockColumns][];
        int blockIndex = 0;
        for (int iBlock = 0; iBlock < blockRows; ++iBlock) {
            final int pStart = iBlock * BLOCK_SIZE;
            final int pEnd = Math.min(pStart + BLOCK_SIZE, rows);
            final int iHeight = pEnd - pStart;
            for (int jBlock = 0; jBlock < blockColumns; ++jBlock) {
                final int qStart = jBlock * BLOCK_SIZE;
                final int qEnd = Math.min(qStart + BLOCK_SIZE, columns);
                final int jWidth = qEnd - qStart;
                blocks[blockIndex] = new double[iHeight * jWidth];
                ++blockIndex;
            }
        }
        return blocks;
    }

    /**
     * Get the inverse of the decomposed matrix.
     *
     * @return the inverse matrix.
     * @throws SingularMatrixException if the decomposed matrix is singular.
     */
    public double[][] getInverse() {
        return solve(ArraysUtils.eye(qrt[0].length));
    }
}