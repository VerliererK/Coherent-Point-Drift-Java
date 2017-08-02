package com.cklin.cpd;

/**
 * The probability matrix for the CPD algorithm
 */
class ProbabilityMatrix {
    /**
     * M x D matrix of assignment
     */
    double[][] PX;
    /**
     * M x 1 vector, sum(P1ij,j=1..N)
     */
    double[] P1;
    /**
     * N x 1 vector, sum(P1ij',i=1..M);
     */
    double[] Pt1;
    /**
     * Summation over all elements of PX
     */
    double E;

    ProbabilityMatrix(double[][] X, double[][] Y) {
        int N = X.length;
        int D = X[0].length;
        int M = Y.length;
        PX = new double[M][D];
        P1 = new double[M];
        Pt1 = new double[N];
        E = 0;
    }
}
