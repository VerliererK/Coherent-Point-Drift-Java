package com.cklin.math;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ArraysUtils {
    //Create

    /**
     * parse string to double[][] array(must be Two-dimensional array also called Matrix), string syntax like
     * <p>"1 2 3; 4 5 6; 7 8 9;"  separate columns by ";" separate rows by space or ","
     *
     * @param a string
     * @return double[][] array
     */
    public static double[][] create(String a) {
        a = a.replaceAll(",", "");
        a = a.trim().replaceAll(" +", " ");
        String[] rowString = a.split(";");
        double[][] A = new double[rowString.length][];

        Scanner scanner = new Scanner("");
        for (int i = 0; i < A.length; i++) {
            scanner = new Scanner(rowString[i]);
            List<Double> scanResult = new ArrayList<Double>();
            while (scanner.hasNextDouble())
                scanResult.add(scanner.nextDouble());
            A[i] = toPrimitive(scanResult);
            if (i > 0 && A[i].length != A[i - 1].length) {
                System.err.println("Dimensions of matrices being concatenated are not consistent.");
                scanner.close();
                return null;
            }
        }
        scanner.close();
        return A;
    }

    public static double[][] zeros(int rows, int columns) {
        return new double[rows][columns];
    }

    public static double[][] zeros(int n) {
        return zeros(n, n);
    }

    public static double[][] ones(int rows, int columns) {
        return ones(rows, columns, 1);
    }

    public static double[][] ones(int n) {
        return ones(n, n);
    }

    public static double[][] ones(int rows, int columns, double scalar) {
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                result[i][j] = scalar;
        return result;
    }

    public static double[][] rand(int n) {
        return rand(n, n);
    }

    public static double[][] rand(int rows, int columns) {
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                result[i][j] = Math.random();
        return result;
    }

    public static double[][] eye(int n) {
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++)
            result[i][i] = 1;
        return result;
    }

    public static double[][] magic(int n) {
        if (n % 4 == 0)
            return magic_4n(n);
        else if (n % 2 == 0)
            return magic_even(n);
        else
            return magic_odd(n);
    }

    private static double[][] magic_4n(int n) {
        double[][] square = new double[n + 1][n + 1];

        for (int j = 1; j <= n; j++) {
            for (int i = 1; i <= n; i++) {
                if (j % 4 == i % 4 || (j % 4 + i) % 4 == 1)
                    square[i][j] = (n + 1 - i) * n - j + 1;
                else
                    square[i][j] = (i - 1) * n + j;
            }
        }

        double[][] matrix = new double[n][n];

        for (int k = 0; k < matrix.length; k++) {
            for (int l = 0; l < matrix[0].length; l++) {
                matrix[k][l] = square[k + 1][l + 1];
            }
        }
        return matrix;
    }

    private static double[][] magic_even(int N) {
        double[][] square = new double[N][N];
        int n = N / 2;
        int row = 0;
        int column = n / 2;
        for (int count = 1; count <= n * n; count++) {
            square[row][column] = count;            // ??A
            square[row + n][column + n] = count + n * n;  // ??B
            square[row][column + n] = count + 2 * n * n;  // ??C
            square[row + n][column] = count + 3 * n * n;  // ??D
            if (count % n == 0)
                row++;
            else {
                row = (row == 0) ? n - 1 : row - 1;
                column = (column == n - 1) ? 0 : column + 1;
            }
        }
        //exchange
        n = N;
        int m = n / 4;
        int m1 = m - 1;

        for (int i = 0; i < n / 2; i++) {
            if (i != m) {
                for (int j = 0; j < m; j++)          //處理規則1
                    magic_swap(square, i, j, n / 2 + i, j);
                for (int j = 0; j < m1; j++)         //處理規則2
                    magic_swap(square, i, n - 1 - j, n / 2 + i, n - 1 - j);
            } else {  //處理規則3
                for (int j = 1; j <= m; j++)
                    magic_swap(square, m, j, n / 2 + m, j);
                for (int j = 0; j < m1; j++)
                    magic_swap(square, m, n - 1 - j, n / 2 + m, n - 1 - j);
            }
        }
        return square;
    }

    private static void magic_swap(double[][] number, int i, int j, int k, int l) {
        double t = number[i][j];
        number[i][j] = number[k][l];
        number[k][l] = t;
    }

    private static double[][] magic_odd(int n) {
        double[][] square = new double[n + 1][n + 1];

        for (int i = 0, j = (n + 1) / 2, key = 1; key <= n * n; key++) {
            if ((key % n) == 1) i++;
            else {
                i--;
                j++;
            }

            if (i == 0) i = n;
            if (j > n) j = 1;

            square[i][j] = key;
        }

        double[][] matrix = new double[n][n];

        for (int k = 0; k < matrix.length; k++) {
            for (int l = 0; l < matrix[0].length; l++) {
                matrix[k][l] = square[k + 1][l + 1];
            }
        }

        return matrix;
    }

    public static double[] diag(double[][] A) {
        int count = Math.min(rowCount(A), colCount(A));
        double[] result = new double[count];
        for (int i = 0; i < count; i++)
            result[i] = A[i][i];
        return result;
    }

    public static double[][] diag(double[] A) {
        int count = A.length;
        double[][] result = zeros(count);
        for (int i = 0; i < count; i++)
            result[i][i] = A[i];
        return result;
    }

    public static double[][] blkdiag(double[][]... ds) {
        int row = 0, col = 0;
        for (double[][] d : ds) {
            row += d.length;
            col += d[0].length;
        }
        double[][] result = new double[row][col];

        int k = 0, l = 0;
        for (double[][] d : ds) {
            for (int i = 0; i < d.length; i++)
                for (int j = 0; j < d[0].length; j++)
                    result[k + i][l + j] = d[i][j];
            k += d.length;
            l += d[0].length;
        }
        return result;
    }

    //Deep Copy
    //1d array of value type can use clone() to deep copy
    public static double[][] copy(double[][] A) {
        double[][] B = new double[A.length][];
        for (int i = 0; i < B.length; i++)
            B[i] = A[i].clone();
        return B;
    }

    public static int[][] copy(int[][] A) {
        int[][] B = new int[A.length][];
        for (int i = 0; i < B.length; i++)
            B[i] = A[i].clone();
        return B;
    }

    //Create Grid
    public static double[] linspace(double x1, double x2, int n) {
        double[] result = new double[n];
        for (int i = 0; i < n; i++)
            result[i] = x1 + (double) i * (x2 - x1) / ((double) n - 1);
        return result;
    }

    public static int[] arithmetic(int start, int diff, int end) {
        int length = (int) (Math.floor((end - start) / diff) + 1);
        if (length <= 0)
            return null;
        int[] A = new int[length];
        for (int i = 0; i < length; i++)
            A[i] = start + diff * i;
        return A;
    }

    public static double[] arithmetic(double start, double diff, double end) {
        int length = (int) (Math.floor((end - start) / diff) + 1);
        double[] A = new double[length];
        for (int i = 0; i < length; i++)
            A[i] = start + diff * i;
        return A;
    }

    //Combine

    /**
     * @return if dimension==1 [d, d, ...]
     * if dimension==2 [d; d; ...]
     */
    public static double[][] cat(int dimension, double[][]... ds) {
        if (dimension == 1) {
            int rows = ds[0].length;
            int cols = 0;
            for (int i = 0; i < ds.length; i++)
                cols += ds[i][0].length;
            double[][] result = new double[rows][cols];
            int k = 0;
            for (int i = 0; i < ds.length; i++) {
                for (int j = 0; j < ds[0].length; j++)
                    System.arraycopy(ds[i][j], 0, result[j], k, ds[i][0].length);
                k += ds[i][0].length;
            }
            return result;
        } else if (dimension == 2) {
            int rows = 0;
            int cols = ds[0][0].length;
            for (int i = 0; i < ds.length; i++)
                rows += ds[i].length;
            double[][] result = new double[rows][cols];

            int k = 0;

            for (int i = 0; i < ds.length; i++)
                for (int j = 0; j < ds[i].length; j++)
                    System.arraycopy(ds[i][j], 0, result[k++], 0, ds[i][j].length);
            return result;
        } else
            return null;
    }

    public static int[][] cat(int dimension, int[][]... ds) {
        if (dimension == 1) {
            int rows = ds[0].length;
            int cols = 0;
            for (int i = 0; i < ds.length; i++)
                cols += ds[i][0].length;
            int[][] result = new int[rows][cols];
            int k = 0;
            for (int i = 0; i < ds.length; i++) {
                for (int j = 0; j < ds[0].length; j++)
                    System.arraycopy(ds[i][j], 0, result[j], k, ds[i][0].length);
                k += ds[i][0].length;
            }
            return result;
        } else if (dimension == 2) {
            int rows = 0;
            int cols = ds[0][0].length;
            for (int i = 0; i < ds.length; i++)
                rows += ds[i].length;
            int[][] result = new int[rows][cols];

            int k = 0;

            for (int i = 0; i < ds.length; i++)
                for (int j = 0; j < ds[i].length; j++)
                    System.arraycopy(ds[i][j], 0, result[k++], 0, ds[i][j].length);
            return result;
        } else
            return null;
    }

    public static int[] cat(int[]... ds) {
        int length = 0;
        for (int i = 0; i < ds.length; i++)
            if (ds[i] != null)
                length += ds[i].length;

        int[] result = new int[length];
        int k = 0;
        for (int i = 0; i < ds.length; i++) {
            if (ds[i] != null) {
                System.arraycopy(ds[i], 0, result, k, ds[i].length);
                k += ds[i].length;
            }
        }
        return result;
    }

    //Size
    public static int rowCount(double[][] A) {
        return A.length;
    }

    public static int colCount(double[][] A) {
        return A[0].length;
    }

    public static int length(double[][] A) {
        return (rowCount(A) > colCount(A)) ? rowCount(A) : colCount(A);
    }

    public static int[] size(double[][] A) {
        return new int[]{rowCount(A), colCount(A)};
    }

    public static int numel(double[][] A) {
        return rowCount(A) * colCount(A);
    }

    //Reshape and Rearrange
    public static double[] getRow(double[][] A, int n) {
        return A[n].clone();
    }

    /**
     * @param A
     * @param n : 0 ~ length-1
     * @return A(n, :)
     */
    public static double[][] getRows(double[][] A, int... n) {
        double[][] result = new double[n.length][A[0].length];
        for (int i = 0; i < n.length; i++)
            result[i] = A[n[i]].clone();
        return result;
    }

    /**
     * @param A
     * @param n : 0 ~ length-1
     * @return A(:, n)
     */
    public static double[] getCol(double[][] A, int n) {
        double[] result = new double[A.length];
        for (int i = 0; i < A.length; i++)
            result[i] = A[i][n];
        return result;
    }

    public static double[][] getCols(double[][] A, int... n) {
        return transpose(getRows(transpose(A), n));
    }

    /**
     * @param A
     * @param indexs
     * @return A(indexes, :) = []
     */
    public static double[][] removeRows(double[][] A, int... indexes) {
        double[][] B = new double[A.length - indexes.length][A[0].length];
        int k = 0;
        outerloop:
        for (int i = 0; i < A.length; i++) {
            for (int index : indexes)
                if (i == index)
                    continue outerloop;
            B[k] = A[i].clone();
            k++;
        }
        return B;
    }

    public static double[][] transpose(double[][] A) {
        double[][] result = new double[A[0].length][A.length];
        for (int i = 0; i < result.length; i++)
            for (int j = 0; j < result[0].length; j++)
                result[i][j] = A[j][i];
        return result;
    }

    public static double[] to1DArray(double[][] A) {
        double[] result = new double[A.length * A[0].length];
        int row = 0, col = 0;
        for (double[] aa : A) {
            col = 0;
            for (double a : aa)
                result[(col++) * A.length + row] = a;
            row++;
        }
        return result;
    }

    public static double[][] reshape(double[] A, int m, int n) {
        double[][] result = new double[m][n];
        int k = 0;
        for (int i = 0; i < result[0].length; i++)
            for (int j = 0; j < result.length; j++)
                result[j][i] = A[k++];
        return result;
    }

    public static double[][] reshape(double[][] A, int m, int n) {
        double[] a = to1DArray(A);
        double[][] result = new double[m][n];
        int k = 0;
        for (int i = 0; i < result[0].length; i++)
            for (int j = 0; j < result.length; j++)
                result[j][i] = a[k++];
        return result;
    }

    public static int[] flip(int[] A) {
        int[] result = A.clone();
        for (int i = 0; i < result.length / 2; i++) {
            int temp = result[i];
            result[i] = result[result.length - 1 - i];
            result[result.length - 1 - i] = temp;
        }
        return result;
    }

    public static double[] flip(double[] A) {
        double[] result = A.clone();
        for (int i = 0; i < result.length / 2; i++) {
            double temp = result[i];
            result[i] = result[result.length - 1 - i];
            result[result.length - 1 - i] = temp;
        }
        return result;
    }

    public static double[][] flip(double[][] A) {
        double[][] result = new double[A.length][A[0].length];
        int k = A.length - 1;
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++)
                result[k][j] = A[i][j];
            k--;
        }
        return result;
    }

    public static double[][] fliplr(double[][] A) {
        double[][] result = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++) {
            int k = A[0].length - 1;
            for (int j = 0; j < A[0].length; j++)
                result[i][k--] = A[i][j];
        }
        return result;
    }

    //Indexing
    public static int[] find(int[] A, int value) {
        ArrayList<Integer> Index = new ArrayList<Integer>();
        int k = 0;
        for (int a : A) {
            if (a == value)
                Index.add(k);
            k++;
        }
        int[] index = new int[Index.size()];
        for (int i = 0; i < index.length; i++)
            index[i] = Index.get(i);
        return index;
    }

    public static int[] find(double[] A, double value) {
        ArrayList<Integer> index = new ArrayList<Integer>();
        int k = 0;
        for (double a : A) {
            if (a == value)
                index.add(k);
            k++;
        }
        return toPrimitive(index.toArray(new Integer[0]));
    }

    public static int[] find(double[][] A, double value) {
        ArrayList<Integer> index = new ArrayList<Integer>();
        int k = 0;
        for (double[] aa : A)
            for (double a : aa) {
                if (a == value)
                    index.add(k);
                k++;
            }
        return toPrimitive((Integer[]) index.toArray());
    }

    public static int[] toPrimitive(Integer[] array) {
        int[] result = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].intValue();
        }
        return result;
    }

    public static double[] toPrimitive(Double[] array) {
        double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].intValue();
        }
        return result;
    }

    public static double[] toPrimitive(List<Double> list) {
        double[] result = new double[list.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = list.get(i);
        }
        return result;
    }

    public static int[] sub2ind(int matrixRowSize, int matrixColSize, final double[] rowSub, final double[] colSub) {
        if (rowSub.length != colSub.length)
            throw new IllegalArgumentException("The subscript vectors must all be of the same size");
        int[] linearInd = new int[rowSub.length];
        for (int i = 0; i < linearInd.length; i++)
            linearInd[i] = sub2ind(matrixRowSize, matrixColSize, (int) Math.round(rowSub[i]), (int) Math.round(colSub[i]));
        return linearInd;
    }

    public static int[] sub2ind(int matrixRowSize, int matrixColSize, int[] rowSub, int[] colSub) {
        if (rowSub.length != colSub.length)
            throw new IllegalArgumentException("The subscript vectors must all be of the same size");
        int[] linearInd = new int[rowSub.length];
        for (int i = 0; i < linearInd.length; i++)
            linearInd[i] = sub2ind(matrixRowSize, matrixColSize, rowSub[i], colSub[i]);
        return linearInd;
    }

    public static int sub2ind(int matrixRowSize, int matrixColSize, int rowSub, int colSub) {
        if (rowSub < matrixRowSize && colSub < matrixColSize)
            return --colSub * matrixRowSize + --rowSub;
        else
            throw new IllegalArgumentException("Out of range subscript");
    }

    //Print
    public static DecimalFormat printFormat = new DecimalFormat("0.0000");

    public static void print(int[] A) {
        if (A == null) {
            System.out.println(A);
            return;
        }
        String result = "Int Matrix: " + "1 x " + A.length + "\n";
        for (int i = 0; i < A.length; i++) {
            result += A[i];
            if (i != A.length - 1)
                result += ", ";
        }
        System.out.println(result);
    }

    public static void print(int[][] A) {
        if (A == null) {
            System.out.println(A);
            return;
        }
        String result = "Int Matrix: " + A.length + " x " + A[0].length;
        for (int i = 0; i < A.length; i++) {
            result += "\n  ";
            for (int j = 0; j < A[0].length; j++) {
                result += A[i][j];
                if (j != A[0].length - 1)
                    result += ", ";
            }
        }
        System.out.println(result);
    }

    public static void print(double[] A) {
        if (A == null) {
            System.out.println(A);
            return;
        }
        String result = "Double Matrix: " + "1 x " + A.length + "\n";
        for (int i = 0; i < A.length; i++) {
            result += printFormat.format(A[i]);
            if (i != A.length - 1)
                result += ", ";
        }
        System.out.println(result);
    }

    public static void print(double[][] A) {
        if (A == null) {
            System.out.println(A);
            return;
        }
        String result = "Double Matrix: " + A.length + " x " + A[0].length;
        for (int i = 0; i < A.length; i++) {
            result += "\n  ";
            for (int j = 0; j < A[0].length; j++) {
                result += printFormat.format(A[i][j]);
                if (j != A[0].length - 1)
                    result += ", ";
            }
        }
        System.out.println(result);
    }
}
