package com.cklin.math;

import java.util.stream.IntStream;

public class ArraysMath {
    public final static int row = 0;
    public final static int col = 1;

    //Elementary Math
    public static double[] add(double[] A, double[] B) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++)
            C[i] = A[i] + B[i];
        return C;
    }

    public static double[][] add(double[][] A, double scalar) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] + scalar;
        return C;
    }

    public static double[][] add(double[][] A, double[] a, int dim) {
        double[][] C = new double[A.length][A[0].length];
        if (dim == row)
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] + a[j];
        else
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] + a[i];
        return C;
    }

    public static double[][] add(double[][] A, double[][] B) {
        checkSameSize(A, B);
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    public static double[] sub(double[] A, double[] B) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++)
            C[i] = A[i] - B[i];
        return C;
    }

    public static double[][] sub(double[][] A, double scalar) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] - scalar;
        return C;
    }

    public static double[][] sub(double[][] A, double[] a, int dim) {
        double[][] C = new double[A.length][A[0].length];
        if (dim == row)
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] - a[j];
        else
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] - a[i];
        return C;
    }

    public static double[][] sub(double[][] A, double[][] B) {
        checkSameSize(A, B);
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    public static double[] mul(double[] A, double[] B) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++)
            C[i] = A[i] * B[i];
        return C;
    }

    public static double[][] mul(double[][] A, double scalar) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] * scalar;
        return C;
    }

    public static double[][] mul(double[][] A, double[] a, int dim) {
        double[][] C = new double[A.length][A[0].length];
        if (dim == row)
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] * a[j];
        else
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] * a[i];
        return C;
    }

    /**
     * A .* B Element-wise multiplication
     *
     * @param A
     * @param B
     * @return
     */
    public static double[][] mul(double[][] A, double[][] B) {
        checkSameSize(A, B);
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] * B[i][j];
        return C;
    }

    /**
     * Matrix Multiplication
     *
     * @param A m &times n
     * @param B n &times p
     * @return double[][] m &times p
     */
    public static double[][] mmul(double[][] A, double[][] B) {
        checkCanMultiply(A, B);
        final int m = A.length;
        final int p = B[0].length;
        final int n = A[0].length;

        final double[][] outData = new double[m][p];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                double temp = A[i][j];
                for (int k = 0; k < p; k++) {
                    outData[i][k] += temp * B[j][k];
                }
            }
        }
        return outData;
    }

    /**
     * Parallel Matrix Multiplication
     *
     * @param A m &times n
     * @param B n &times p
     * @return double[][] m &times p
     */
    public static double[][] mmul_P(double[][] A, double[][] B) {
        checkCanMultiply(A, B);
        final int m = A.length;
        final int p = B[0].length;
        final int n = A[0].length;

        final double[][] outData = new double[m][p];
        IntStream.range(0, m).parallel().forEach(i -> {
            for (int j = 0; j < n; j++) {
                double temp = A[i][j];
                for (int k = 0; k < p; k++) {
                    outData[i][k] += temp * B[j][k];
                }
            }
        });
        return outData;
    }

    public static double[] div(double[] A, double[] B) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++)
            C[i] = A[i] / B[i];
        return C;
    }

    /**
     * A ./ scalar Element-wise division
     *
     * @param A
     * @param scalar
     * @return
     */
    public static double[][] div(double[][] A, double scalar) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] / scalar;
        return C;
    }

    /**
     * scalar ./ A Element-wise division
     *
     * @param scalar
     * @param A
     * @return
     */
    public static double[][] div(double scalar, double[][] A) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = scalar / A[i][j];
        return C;
    }

    /**
     * A ./ a Element-wise division
     *
     * @param A
     * @param a
     * @param dim dimension to operate
     * @return
     */
    public static double[][] div(double[][] A, double[] a, int dim) {
        double[][] C = new double[A.length][A[0].length];
        if (dim == row)
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] / a[j];
        else
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = A[i][j] / a[i];
        return C;
    }

    /**
     * a ./ A Element-wise division
     *
     * @param a
     * @param A
     * @param dim dimension to operate
     * @return
     */
    public static double[][] div(double[] a, double[][] A, int dim) {
        double[][] C = new double[A.length][A[0].length];
        if (dim == row)
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = a[j] / A[i][j];
        else
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    C[i][j] = a[i] / A[i][j];
        return C;
    }

    /**
     * A ./ B Element-wise division
     *
     * @param A
     * @param B
     * @return
     */
    public static double[][] div(double[][] A, double[][] B) {
        checkSameSize(A, B);
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = A[i][j] / B[i][j];
        return C;
    }

    public static double[][] neg(double[][] A) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = -A[i][j];
        return C;
    }

    public static double[][] abs(double[][] A) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = Math.abs(A[i][j]);
        return C;
    }

    public static double[][] round(double[][] A) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = Math.round(A[i][j]);
        return C;
    }

    /**
     * A.^scalar Element-wise power
     *
     * @param A
     * @param scalar
     * @return
     */
    public static double[][] pow(double[][] A, double scalar) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = Math.pow(A[i][j], scalar);
        return C;
    }

    public static double[] sign(double[] A) {
        double[] C = new double[A.length];
        for (int i = 0; i < C.length; i++)
            C[i] = Math.signum(A[i]);
        return C;
    }

    public static double[][] sign(double[][] A) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++)
                C[i][j] = Math.signum(A[i][j]);
        return C;
    }

    public static int[] diff(int[] A) {
        int[] result = new int[A.length - 1];
        for (int i = 0; i < result.length; i++)
            result[i] = A[i + 1] - A[i];
        return result;
    }

    public static double[] diff(double[] A) {
        double[] result = new double[A.length - 1];
        for (int i = 0; i < result.length; i++)
            result[i] = A[i + 1] - A[i];
        return result;
    }

    public static boolean equal(double[] A, double[] B, double threhold) {
        boolean eq = false;
        for (int i = 0; i < A.length; i++) {
            eq = Math.abs(A[i] - B[i]) <= threhold;
            if (eq == false) return eq;
        }
        return eq;
    }

    public static boolean equal(double[][] A, double[][] B, double threhold) {
        boolean eq = false;
        for (int i = 0; i < A.length; i++)
            for (int j = 0; j < A[0].length; j++) {
                eq = Math.abs(A[i][j] - B[i][j]) <= threhold;
                if (eq == false) return eq;
            }
        return eq;
    }

    //Check Matrix Dimension
    private static void checkSameSize(double[][] A, double[][] B) {
        if (A.length != B.length || A[0].length != B[0].length)
            throw new IllegalArgumentException("A size must same as B");
    }

    private static void checkCanMultiply(double[][] A, double[][] B) {
        if (A[0].length != B.length)
            throw new IllegalArgumentException("illegal matrix dimensions");
    }

    //Statistics
    public static int sum(int[] A) {
        int sum = 0;
        for (int a : A)
            sum += a;
        return sum;
    }

    public static double sum(double[] A) {
        double sum = 0;
        for (double a : A)
            sum += a;
        return sum;
    }

    public static double[] sum(double[][] A) {
        if (A.length == 1)
            return sum(A, row);
        else
            return sum(A, col);
    }

    public static double[] sum(double[][] A, int dim) {
        if (dim == col) {
            double[] sum = new double[A[0].length];
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    sum[j] += A[i][j];
            return sum;
        } else if (dim == row) {
            double[] sum = new double[A.length];
            for (int i = 0; i < A.length; i++)
                for (int j = 0; j < A[0].length; j++)
                    sum[i] += A[i][j];
            return sum;
        } else
            return null;
    }

    public static double mean(int[] A) {
        return (double) sum(A) / A.length;
    }

    public static double mean(double[] A) {
        return sum(A) / A.length;
    }

    public static double[] mean(double[][] A) {
        double[] result;
        if (A.length == 1)
            result = mean(A, row);
        else
            result = mean(A, col);
        return result;
    }

    public static double[] mean(double[][] A, int dim) {
        double[] result = sum(A, dim);
        int N = (result.length == A.length) ? A[0].length : A.length;
        for (int j = 0; j < result.length; j++)
            result[j] /= N;
        return result;
    }

    public static double var(double[] A) {
        double var = Double.NaN;
        double mean = mean(A);
        double length = A.length;
        if (length == 1) {
            var = 0.0;
        } else if (length > 1) {
            double accum = 0.0;
            double dev = 0.0;
            double accum2 = 0.0;
            for (int i = 0; i < length; i++) {
                dev = A[i] - mean;
                accum += dev * dev;
                accum2 += dev;
            }
            var = (accum - (accum2 * accum2 / length)) / (length - 1.0);
        }
        return var;
    }

    public static double std(double[] A) {
        return Math.sqrt(var(A));
    }

    public static double max(double[] A) {
        int index = maxIndex(A);
        return A[index];
    }

    public static double[] max(double[][] A) {
        int[] index = maxIndex(A);
        double[] max = new double[index.length];
        for (int i = 0; i < index.length; i++)
            max[i] = A[i][index[i]];
        return max;
    }

    public static int maxIndex(double[] A) {
        double max = A[0];
        int index = 0;
        for (int i = 0; i < A.length; i++)
            if (A[i] > max) {
                index = i;
                max = A[i];
            }
        return index;
    }

    public static int[] maxIndex(double[][] A) {
        int[] index = new int[A.length];
        for (int i = 0; i < A.length; i++)
            index[i] = maxIndex(A[i]);
        return index;
    }

    public static double min(double[] A) {
        int index = minIndex(A);
        return A[index];
    }

    public static double[] min(double[][] A) {
        int[] index = minIndex(A);
        double[] max = new double[index.length];
        for (int i = 0; i < index.length; i++)
            max[i] = A[i][index[i]];
        return max;
    }

    public static int minIndex(double[] A) {
        double min = A[0];
        int index = 0;
        for (int i = 0; i < A.length; i++)
            if (A[i] < min) {
                index = i;
                min = A[i];
            }
        return index;
    }

    public static int[] minIndex(double[][] A) {
        int[] index = new int[A.length];
        for (int i = 0; i < A.length; i++)
            index[i] = minIndex(A[i]);
        return index;
    }

    //Linear Algebra
    public static double dot(double[] A, double[] B) {
        if (A.length != B.length)
            throw new IllegalArgumentException("A size must same as B");
        double C = 0;
        for (int i = 0; i < A.length; i++)
            C += A[i] * B[i];
        return C;
    }

    public static double[] cross(double[] A, double[] B) {
        if (A.length != 3)
            throw new IllegalArgumentException("A must be length of 3");
        if (A.length != B.length)
            throw new IllegalArgumentException("A size must same as B");
        double[] C = new double[3];
        C[0] = A[1] * B[2] - A[2] * B[1];
        C[1] = -(A[0] * B[2] - A[2] * B[0]);
        C[2] = A[0] * B[1] - A[1] * B[0];
        return C;
    }

    public static double norm(double[] A) {
        double C = 0;
        for (int i = 0; i < A.length; i++)
            C += A[i] * A[i];
        C = Math.sqrt(C);
        return C;
    }

    public static double[] norm(double[][] A) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++)
            C[i] = norm(A[i]);
        return C;
    }

    public static double[] normr(double[] A) {
        double[] C = new double[A.length];
        double norm = norm(A);
        for (int i = 0; i < C.length; i++)
            C[i] = A[i] / norm;
        return C;
    }

    public static double[][] normr(double[][] A) {
        double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++)
            C[i] = normr(A[i]);
        return C;
    }

    public static double trace(double[][] A) {
        return sum(ArraysUtils.diag(A));
    }

    public static double det(double[][] A) {
        LUDecomposition lu = new LUDecomposition(A);
        return lu.getDeterminant();
    }

    public static double[][] inv(double[][] A) {
        LUDecomposition lu = new LUDecomposition(A);
        return lu.getInverse();
    }

    /**
     * Solve systems of linear equations Ax = B for x
     * <p>
     * same as mldivide A\B
     *
     * @param A
     * @param B
     * @return x
     */
    public static double[][] solve(double[][] A, double[][] B) {
        if (A.length == A[0].length)    //Square Matrix
            return new LUDecomposition(A).solve(B);
        else
            return new QRDecomposition(A).solve(B);
    }

    /**
     * Solve systems of linear equations xA = B for x
     * <p>
     * same as mrdivide B/A
     *
     * @param A
     * @param B
     * @return x
     */
    public static double[][] solveTranspose(double[][] A, double[][] B) {
        return ArraysUtils.transpose(solve(ArraysUtils.transpose(A), ArraysUtils.transpose(B)));
    }

    /**
     * A\B, Solve systems of linear equations Ax = B for x
     *
     * @param A
     * @param B
     * @return x
     */
    public static double[][] mldivide(double[][] A, double[][] B) {
        return solve(A, B);
    }

    /**
     * B/A, Solve systems of linear equations xA = B for x
     *
     * @param A
     * @param B
     * @return x
     */
    public static double[][] mrdivide(double[][] A, double[][] B) {
        return solveTranspose(A, B);
    }

    //Polynomial
    //TODO Polynomial

    //Others
    public static double distance(double[] A, double[] B) {
        return norm(sub(A, B));
    }

    public static double[] smooth(double[] A) {
        double[] B = new double[A.length];
        B[0] = A[0];
        B[1] = (A[0] + A[1] + A[2]) / 3;
        for (int i = 2; i < A.length - 2; i++)
            B[i] = (A[i - 2] + A[i - 1] + A[i] + A[i + 1] + A[i + 2]) / 5;
        B[A.length - 2] = (A[A.length - 3] + A[A.length - 2] + A[A.length - 1]) / 3;
        B[A.length - 1] = A[A.length - 1];
        return B;
    }

    public static int[] findpeak(double[] A, String type) {
        double[] B = diff(A);
        double[] C = sign(B);
        double[] D = diff(C);
        int[] locations;
        switch (type) {
            case "max":
                locations = ArraysUtils.find(D, -2);
                for (int i = 0; i < locations.length; i++)
                    locations[i]++;
                return locations;
            case "min":
                locations = ArraysUtils.find(D, 2);
                for (int i = 0; i < locations.length; i++)
                    locations[i]++;
                return locations;
            default:
                throw new IllegalArgumentException("What kind of peak type (max or min) do you want to search?");
        }
    }

    public static double[][] interp2(double[][] src, int fx, int fy) {
        double[][] out = new double[fy * (src.length - 1) + 1][fx * (src[0].length - 1) + 1];
        for (int i = 0; i < out.length; i++)
            for (int j = 0; j < out[0].length; j++) {
                int p1_x = i / fy;
                int p1_y = j / fx;

                int p2_x = p1_x + 1;
                int p2_y = p1_y;

                int p3_x = p1_x;
                int p3_y = p1_y + 1;

                int p4_x = p1_x + 1;
                int p4_y = p1_y + 1;

                if (p2_x >= src.length)
                    p2_x--;
                if (p4_x >= src.length)
                    p4_x--;
                if (p3_y >= src[0].length)
                    p3_y--;
                if (p4_y >= src[0].length)
                    p4_y--;

                double alpha = (double) i / (double) fy - (double) p1_x;
                double beta = (double) j / (double) fx - (double) p1_y;

                double value =
                        (1 - alpha) * (1 - beta) * src[p1_x][p1_y] +
                                alpha * (1 - beta) * src[p2_x][p2_y] +
                                (1 - alpha) * beta * src[p3_x][p3_y] +
                                alpha * beta * src[p4_x][p4_y];
                out[i][j] = value;
            }
        return out;
    }

    //Coordinate
    public static double[][] zyx2r(double A1, double A2, double A3) {
        return zyx2r(new double[]{A1, A2, A3});
    }

    public static double[][] zyx2r(double A[]) {
        for (int i = 0; i < A.length; i++)
            A[i] = A[i] * Math.PI / 180;
        double R[][] = new double[3][3];
        R[0][0] = Math.cos(A[2]) * Math.cos(A[1]);
        R[0][1] = -Math.sin(A[2]) * Math.cos(A[0]) + Math.cos(A[2]) * Math.sin(A[1]) * Math.sin(A[0]);
        R[0][2] = Math.sin(A[2]) * Math.sin(A[0]) + Math.cos(A[2]) * Math.sin(A[1]) * Math.cos(A[0]);
        R[1][0] = Math.sin(A[2]) * Math.cos(A[1]);
        R[1][1] = Math.cos(A[2]) * Math.cos(A[0]) + Math.sin(A[2]) * Math.sin(A[1]) * Math.sin(A[0]);
        R[1][2] = -Math.cos(A[2]) * Math.sin(A[0]) + Math.sin(A[2]) * Math.sin(A[1]) * Math.cos(A[0]);
        R[2][0] = -Math.sin(A[1]);
        R[2][1] = Math.cos(A[1]) * Math.sin(A[0]);
        R[2][2] = Math.cos(A[1]) * Math.cos(A[0]);
        return R;
    }
}