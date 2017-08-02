package com.cklin.cpd;

import static com.cklin.math.ArraysUtils.colCount;
import static com.cklin.math.ArraysUtils.copy;
import static com.cklin.math.ArraysUtils.rowCount;
import static com.cklin.math.ArraysUtils.transpose;

import com.cklin.math.ArraysMath;
import com.cklin.math.ArraysUtils;

import com.cklin.math.SingularValueDecomposition;
import test.Plot;

public class Rigid {
    private static final double eps = 1e-12;

    public static double[][] rigid(CPD cpd) {
        int max_it = cpd.opt.max_iter;
        double tol = cpd.opt.tol;
        double outliers = cpd.opt.outliers;
        //boolean corresp = cpd.opt.corresp;
        boolean fgt = cpd.opt.fgt;

        int N = rowCount(cpd.X);
        int D = colCount(cpd.X);
        int M = rowCount(cpd.Y);

        // Initialize sigma and Y
        double sigma2 = Utils.getSigma2(cpd);
        //double sigma2_init = sigma2;
        double[][] T = copy(cpd.Y);
        double s = 1;
        double[][] R = ArraysUtils.eye(D);

        //Plot
        Plot plot = null;
        if (cpd.opt.plot) {
            plot = new Plot();
            plot.show();
            plot.addXYSeries("X", cpd.X);
            plot.addXYSeries("T", T);
        }
        // Optimization
        int iter = 0;
        double ntol = tol + 10;
        double L = 0;
        while ((iter < max_it) && (ntol > tol) && (sigma2 > 10 * eps)) {
            double L_old = L;
            // Check wheather we want to use the Fast Gauss Transform
            ProbabilityMatrix P = (fgt) ? Utils.cpd_P(cpd.X, T, sigma2, outliers) :
                    Utils.cpd_P(cpd.X, T, sigma2, outliers);

            L = P.E;
            ntol = Math.abs((L - L_old) / L);
            if (cpd.opt.print)
                System.out.println(
                        "CPD Rigid" + " : dL= " + ntol +
                                ", iter= " + iter + ", sigma2= " + sigma2);
            // Precompute
            double Np = ArraysMath.sum(P.Pt1);
            double[] mu_x = new double[D];
            double[] mu_y = new double[D];
            //mu_x=X'*Pt1/Np (DxN)*(N*1);
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < N; j++)
                    mu_x[i] += cpd.X[j][i] * P.Pt1[j];
                mu_x[i] /= Np;
            }
            // mu_y=Y'*P1/Np (DxM)*(M*1);
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < M; j++)
                    mu_y[i] += cpd.Y[j][i] * P.P1[j];
                mu_y[i] /= Np;
            }

            // Solve for Rotation, scaling, translation and sigma^2
            double[][] A_temp = ArraysMath.mmul(transpose(P.PX), cpd.Y);
            double[][] A = ArraysMath.sub(A_temp, mmul1D(mu_x, mu_y, Np));
            SingularValueDecomposition svd = new SingularValueDecomposition(A);
            double[][] U = svd.getU();
            double[][] S = svd.getS();
            double[][] V = svd.getV();
            double[][] C = ArraysUtils.eye(D);
            // check if we need strictly rotation (no reflections)
            if (cpd.opt.rot)
                C[D - 1][D - 1] = ArraysMath.det(ArraysMath.mmul(U, V));
            R = ArraysMath.mmul(ArraysMath.mmul(U, C), transpose(V));
            // check if estimating scaling as well, otherwise s=1
            cpd.sigma2 = sigma2;
            if (cpd.opt.scale) {
                double traceSC = ArraysMath.trace(ArraysMath.mmul(S, C));
                double temp_YYP1 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(cpd.Y, cpd.Y), P.P1, ArraysMath.col)));
                double temp_mu_y = Np * ArraysMath.sum(ArraysMath.mul(mu_y, mu_y));
                s = traceSC / (temp_YYP1 - temp_mu_y);
                double temp_XXPt1 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(cpd.X, cpd.X), P.Pt1, ArraysMath.col)));
                double temp_mu_x = Np * ArraysMath.sum(ArraysMath.mul(mu_x, mu_x));
                sigma2 = Math.abs(temp_XXPt1 - temp_mu_x - s * traceSC) / (Np * D);
            } else {
                double traceSC = ArraysMath.trace(ArraysMath.mmul(S, C));
                double temp_XXPt1 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(cpd.X, cpd.X), P.Pt1, ArraysMath.col)));
                double temp_mu_x = Np * ArraysMath.sum(ArraysMath.mul(mu_x, mu_x));
                double temp_YYP1 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(cpd.Y, cpd.Y), P.P1, ArraysMath.col)));
                double temp_mu_y = Np * ArraysMath.sum(ArraysMath.mul(mu_y, mu_y));
                sigma2 = Math.abs(temp_XXPt1 - temp_mu_x + temp_YYP1 - temp_mu_y
                        - 2 * traceSC) / (Np * D);
            }
            // t=mu_x-s*R*mu_y;
            double[] t = new double[mu_x.length];
            for (int i = 0; i < t.length; i++) {
                double temp = 0;
                for (int j = 0; j < R[0].length; j++)
                    temp += R[i][j] * mu_y[j];
                t[i] = mu_x[i] - s * temp;
            }

            // Update the GMM centroids
            //T=s*Y*R'+repmat(t',[M 1]);
            T = ArraysMath.mul(ArraysMath.mmul(cpd.Y, transpose(R)), s);
            T = ArraysMath.add(T, t, ArraysMath.row);
            iter = iter + 1;

            if (cpd.opt.plot) {
                plot.setXYSeries("T", T);
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
        System.out.println("CPD registration succesfully completed.");
        return T;
    }

    public static double[][] affine(CPD cpd) {
        int max_it = cpd.opt.max_iter;
        double tol = cpd.opt.tol;
        double outliers = cpd.opt.outliers;
        //boolean corresp = cpd.opt.corresp;
        boolean fgt = cpd.opt.fgt;

        int N = rowCount(cpd.X);
        int D = colCount(cpd.X);
        int M = rowCount(cpd.Y);

        // Initialize sigma and Y
        double sigma2 = Utils.getSigma2(cpd);
        //double sigma2_init = sigma2;
        double[][] T = copy(cpd.Y);

        //Plot
        Plot plot = null;
        if (cpd.opt.plot) {
            plot = new Plot();
            plot.show();
            plot.addXYSeries("X", cpd.X);
            plot.addXYSeries("T", T);
        }
        // Optimization
        int iter = 0;
        double ntol = tol + 10;
        double L = 1;
        while ((iter < max_it) && (ntol > tol) && (sigma2 > 10 * eps)) {
            double L_old = L;
            // Check wheather we want to use the Fast Gauss Transform
            ProbabilityMatrix P = (fgt) ? Utils.cpd_P(cpd.X, T, sigma2, outliers) :
                    Utils.cpd_P(cpd.X, T, sigma2, outliers);
            L = P.E;
            ntol = Math.abs((L - L_old) / L);
            if (cpd.opt.print)
                System.out.println(
                        "CPD Affine" + " : dL= " + ntol +
                                ", iter= " + iter + ", sigma2= " + sigma2);
            // Precompute
            double Np = ArraysMath.sum(P.P1);
            double[] mu_x = new double[D];
            double[] mu_y = new double[D];
            //mu_x=X'*Pt1/Np (DxN)*(N*1);
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < N; j++)
                    mu_x[i] += cpd.X[j][i] * P.Pt1[j];
                mu_x[i] /= Np;
            }
            // mu_y=Y'*P1/Np (DxM)*(M*1);
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < M; j++)
                    mu_y[i] += cpd.Y[j][i] * P.P1[j];
                mu_y[i] /= Np;
            }
            // Solve for parameters
            double[][] B1_temp = ArraysMath.mmul(transpose(P.PX), cpd.Y);
            double[][] B1 = ArraysMath.sub(B1_temp, mmul1D(mu_x, mu_y, Np));

            double[][] B2_temp = ArraysMath.mmul(transpose(ArraysMath.mul(cpd.Y, P.P1, ArraysMath.col)), cpd.Y);
            double[][] B2 = ArraysMath.sub(B2_temp, mmul1D(mu_y, mu_y, Np));
            double[][] B = ArraysMath.mrdivide(B2, B1);
            //B1=PX'*Y - Np*(mu_x*mu_y');
            //B2=(Y.*repmat(P1,1,D))'*Y - Np*(mu_y*mu_y');
            //B=B1/B2; % B= B1 * inv(B2);

            double[] t = new double[D];
            for (int i = 0; i < D; i++) {
                double temp = 0;
                for (int j = 0; j < D; j++)
                    temp += B[i][j] * mu_y[j];
                t[i] = mu_x[i] - temp;
            }
            //t=mu_x-B*mu_y;

            double temp1 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(cpd.X, cpd.X), P.Pt1, ArraysMath.col)));
            double temp2 = Np * ArraysMath.sum(ArraysMath.mul(mu_x, mu_x));
            double temp3 = -1 * ArraysMath.trace(ArraysMath.mmul(B1, transpose(B)));

            cpd.sigma2 = sigma2;
            sigma2 = Math.abs((temp1 + temp2 + temp3) / (Np * D));
            //sigma2=abs(sum(sum(X.^2.*repmat(Pt1,1,D)))- Np*(mu_x'*mu_x) -trace(B1*B'))/(Np*D);

            // abs here to prevent roundoff errors that leads to negative sigma^2 in
            // rear cases

            // Update centroids positioins

            // T=Y*B'+repmat(t',[M 1]);
            double[][] T_temp = ArraysMath.mmul(cpd.Y, transpose(B));
            for (int i = 0; i < T.length; i++)
                for (int j = 0; j < T[0].length; j++)
                    T[i][j] = T_temp[i][j] + t[j];

            iter = iter + 1;

            if (cpd.opt.plot) {
                plot.setXYSeries("T", T);
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
        System.out.println("CPD registration succesfully completed.");
        return T;
    }

    /**
     * Np*(mu_x*mu_y')
     *
     * @param mu_x Dx1
     * @param mu_y 1xD
     * @param Np
     * @return DxD
     */
    private static double[][] mmul1D(double[] mu_x, double[] mu_y, double Np) {
        int N = mu_x.length;
        if (N != mu_y.length)
            throw new IllegalArgumentException("mu_x length must same as mu_y");
        double[][] mu = new double[N][N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                mu[i][j] = mu_x[i] * mu_y[j] * Np;
        return mu;
    }
}
