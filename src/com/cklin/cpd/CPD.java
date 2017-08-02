package com.cklin.cpd;

import com.cklin.math.ArraysMath;
import com.cklin.math.ArraysUtils;

/**
 * align Y onto X.
 */
public class CPD {
    /**
     * registration options
     */
    public Option opt = new Option();
    /**
     * final sigma^2
     */
    public double sigma2;
    /**
     * variables that can be used to re-scale and shift the point
     */
    protected Normal normal;
    /**
     * N x D full 2-D matrices of point-set locations.<br>
     * D is the dimension of point-sets.
     */
    public double[][] X;
    /**
     * M x D full 2-D matrices of point-set locations.<br>
     * D is the dimension of point-sets.
     */
    public double[][] Y;
    /**
     * registered Y point-set
     */
    protected double[][] registeredY;

    /**
     * align Y onto X.
     *
     * @param X N x D full 2-D matrices of point-set locations.
     * @param Y M x D full 2-D matrices of point-set locations.
     */
    public CPD(double[][] X, double[][] Y) {
        this.X = ArraysUtils.copy(X);
        this.Y = ArraysUtils.copy(Y);
    }

    /**
     * @return registered Y point-set
     */
    public double[][] register() {
        if (ArraysUtils.colCount(X) != ArraysUtils.colCount(Y))
            throw new IllegalArgumentException("The dimension of point-sets must be same");
        checkOpt();
        Utils.cpd_register(this);
        return registeredY;
    }

    public double[] interp(double[] point, double[][] originY, double csigma2, double outlier) {
        int N = 1;
        int M = Y.length;
        int D = Y[0].length;
        double ksig = -2.0 * csigma2;
        double outlier_tmp = (outlier * M * Math.pow(-ksig * Math.PI, 0.5 * D)) / ((1 - outlier) * N);

        double[] P = new double[M];
        double[] Dis = new double[M];
        for (int i = 0; i < M; i++) {
            double dis = ArraysMath.distance(point, this.Y[i]);
            Dis[i] = dis * dis;
            P[i] = Math.exp(dis * dis / ksig);
        }
        double Psum = ArraysMath.sum(P) + outlier_tmp;
        for (int i = 0; i < M; i++)
            P[i] /= Psum;
        double[] newPoint = point.clone();
        for (int i = 0; i < M; i++)
            for (int j = 0; j < D; j++) {
                newPoint[j] += P[i] * (registeredY[i][j] - originY[i][j]);
            }
        return newPoint;
    }

    /**
     * @return Correspondence Index
     */
    public int[] getCorrespondence() {
        if (!opt.corresp)
            throw new IllegalArgumentException(
                    "You should set opt.corresp = true, and register again if you modify other options");
        int[] C = Utils.cpd_Pcorrespondence(this);
        return C;
    }

    /**
     * checking for the option
     */
    private void checkOpt() {
        if (this.opt.outliers < 0.0 || this.opt.outliers > 1.0)
            throw new IllegalArgumentException("The weight of noise and outliers must be 0.0~1.0");


        //Non-rigid registration options
        if (this.opt.lambda < 0.0)
            throw new IllegalArgumentException("Regularization weight must greater than 0.0");
        if (this.opt.beta < 0.0)
            throw new IllegalArgumentException("Gaussian smoothing filter size must greater than 0.0");
    }


}
