package com.cklin.cpd;

import static com.cklin.math.ArraysUtils.*;

import com.cklin.math.ArraysMath;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import test.Plot;

class Utils {
    /**
     * Rigid, affine, non-rigid point set registration.
     *
     * @param cpd {@link CPD}
     */
    protected static void cpd_register(CPD cpd) {
        if (cpd.opt.normalize)
            cpd.normal = cpd_normal(cpd);

        double[][] T = null;
        switch (cpd.opt.method) {
            case NonRigid:
                if (cpd.opt.plot)
                    T = cpd_GRBF(cpd);
                else
                    T = NonRigid.GRBF(cpd);
                break;
            case NonRigid_lowRank:
                throw new NotImplementedException();
            case Rigid:
                T = Rigid.rigid(cpd);
                break;
            case Affine:
                T = Rigid.affine(cpd);
                break;
        }

        if (cpd.opt.normalize) {
            for (int i = 0; i < rowCount(T); i++)
                for (int j = 0; j < colCount(T); j++)
                    T[i][j] = T[i][j] * cpd.normal.xScale + cpd.normal.xd[j];
        }
        cpd.registeredY = T;
    }

    /**
     * Normalizes x and y to have zero mean and unit variance.
     *
     * @param cpd {@link CPD}
     * @return {@link Normal}
     */
    private static Normal cpd_normal(CPD cpd) {
        int N = rowCount(cpd.X);
        int M = rowCount(cpd.Y);
        int D = colCount(cpd.X);

        double[] xd = ArraysMath.mean(cpd.X);
        double[] yd = ArraysMath.mean(cpd.Y);


        for (int i = 0; i < N; i++)
            for (int j = 0; j < D; j++)
                cpd.X[i][j] -= xd[j];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < D; j++)
                cpd.Y[i][j] -= yd[j];

        double xScale = 0, yScale = 0;
        for (int i = 0; i < cpd.X.length; i++) {
            for (final double val : cpd.X[i]) {
                xScale += val * val;
            }
        }
        for (int i = 0; i < cpd.Y.length; i++) {
            for (final double val : cpd.Y[i]) {
                yScale += val * val;
            }
        }

        xScale = Math.sqrt(xScale / N);
        yScale = Math.sqrt(yScale / M);

        for (int i = 0; i < N; i++)
            for (int j = 0; j < D; j++)
                cpd.X[i][j] /= xScale;
        for (int i = 0; i < M; i++)
            for (int j = 0; j < D; j++)
                cpd.Y[i][j] /= yScale;

        return new Normal(xd, yd, xScale, yScale);
    }


    /**
     * Construct Gaussian affinity matrix
     *
     * @param X    data points-set
     * @param Y    data points-set
     * @param beta regularization weight
     * @return Gaussian matrix
     */
    protected static double[][] cpd_G(double[][] X, double[][] Y, double beta) {
        double k = -2 * beta * beta;
        int N = rowCount(X);
        int M = rowCount(Y);
        int D = colCount(X);

        double[][] G = new double[N][M];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++) {
                double value = 0;
                for (int d = 0; d < D; d++)
                    value += (X[j][d] - Y[i][d]) * (X[j][d] - Y[i][d]);
                value = Math.exp(value / k);
                G[j][i] = value;
            }
        return G;
    }

    /**
     * Computes the probability matrix for the CPD algorithm
     *
     * @param X       data points-set
     * @param Y       data points-set
     * @param sigma2  standard deviation of Gaussians
     * @param outlier noise weight
     * @return probability matrix
     */
    protected static ProbabilityMatrix cpd_P(double[][] X, double[][] Y, double sigma2, double outlier) {
        int N = rowCount(X);
        int D = colCount(X);
        int M = rowCount(Y);
        int n = 0, m = 0, d = 0;
        double ksig, diff, razn, outlier_tmp, sp;

        ProbabilityMatrix PM = new ProbabilityMatrix(X, Y);

        double[] P = new double[M];
        double[] temp_x = new double[D];
        ;
        ksig = -2.0 * sigma2;
        outlier_tmp = (outlier * M * Math.pow(-ksig * Math.PI, 0.5 * D)) / ((1 - outlier) * N);

        for (n = 0; n < N; n++) {
            sp = 0;
            for (m = 0; m < M; m++) {
                razn = 0;
                for (d = 0; d < D; d++) {
                    diff = X[n][d] - Y[m][d];
                    diff = diff * diff;
                    razn += diff;
                }
                double value = Math.exp(razn / ksig);
                P[m] = value;
                sp += value;
            }
            sp += outlier_tmp;
            PM.Pt1[n] = 1 - outlier_tmp / sp;
            for (d = 0; d < D; d++)
                temp_x[d] = X[n][d] / sp;

            for (m = 0; m < M; m++) {
                PM.P1[m] += P[m] / sp;
                //P1[m, 0] += (P[m, 0] / sp);
                for (d = 0; d < D; d++)
                    PM.PX[m][d] = PM.PX[m][d] + temp_x[d] * P[m];
            }
            PM.E += -Math.log(sp);
        }
        PM.E += D * N * Math.log(sigma2) / 2.0;

        return PM;
    }

    protected static ProbabilityMatrix cpd_P_FGT(double[][] X, double[][] Y,
                                                 double sigma2, double outliers, double sigma2_init) {
        sigma2 = Math.max(0.05, sigma2);
        int N = rowCount(X);
        int D = colCount(X);
        int M = rowCount(Y);

        double hsigma = Math.sqrt(2 * sigma2);
        outliers = (outliers == 0) ? 1e-16 : outliers;

        //FGT parameters
        int e = 9; //Ratio of far field (default e = 10)
        int K = Math.min(N, M);
        K = (int) Math.round(Math.min(K, 50 + sigma2_init / sigma2)); // Number of centers (default K = sqrt(Nx))
        int p = 6; // Order of truncation (default p = 8)

        //computer Pt1 and denomP
        double[] Pt1_w = new double[M];
        for (int i = 0; i < M; i++)
            Pt1_w[i] = 1;
        FGT model = new FGT(transpose(Y),
                Pt1_w, hsigma, e, K, p);
        double[] Kt1 = model.predict(transpose(X), hsigma, e);

        double ndi = outliers / (1 - outliers) * M / N * Math.pow(2 * Math.PI * sigma2, 0.5 * D);
        double[] denomP = new double[Kt1.length];
        for (int i = 0; i < denomP.length; i++)
            denomP[i] = Kt1[i] + ndi;
        double[] Pt1 = new double[denomP.length];
        for (int i = 0; i < Pt1.length; i++)
            Pt1[i] = 1 - ndi / denomP[i];
        //Pt1=Pt1';

        //compute P1
        double[] P1_w = new double[denomP.length];
        for (int i = 0; i < P1_w.length; i++)
            P1_w[i] = 1.0 / denomP[i];
        FGT P1model = new FGT(transpose(X),
                P1_w, hsigma, e, K, p);
        double[] P1 = P1model.predict(transpose(Y), hsigma, e);
        //P1=P1';

        //compute PX
        double[][] PX = new double[D][];
        for (int i = 0; i < D; i++) {
            double[] PX_w = new double[denomP.length];
            for (int j = 0; j < PX_w.length; j++)
                PX_w[j] = X[j][i] / denomP[j];
            //X(:,i)'./denomP
            FGT PXmodel = new FGT(transpose(X),
                    PX_w, hsigma, e, K, p);
            PX[i] = PXmodel.predict(transpose(Y), hsigma, e);
        }
        //PX=PX';
        double denomP_log_sum = 0.0;
        for (int i = 0; i < denomP.length; i++)
            denomP_log_sum += (Math.log(denomP[i]));
        double L = -denomP_log_sum + D * N * Math.log(sigma2) / 2;
        ProbabilityMatrix PM = new ProbabilityMatrix(X, Y);
        PM.E = L;
        PM.PX = transpose(PX);
        PM.Pt1 = Pt1;
        PM.P1 = P1;

        return PM;
    }

    /**
     * @param cpd
     * @return correspondence index
     */
    protected static int[] cpd_Pcorrespondence(CPD cpd) {
        int N = rowCount(cpd.X);
        int D = colCount(cpd.X);
        int M = rowCount(cpd.registeredY);

        int[] C = new int[M];

        int n, m, d;
        double ksig, diff, razn, outlier_tmp, temp_x, sp;
        double[] P = new double[M];
        double[] P1 = new double[M];

        ksig = -2.0 * (cpd.sigma2 + 1e-3);
        outlier_tmp = (cpd.opt.outliers * M * Math.pow(-ksig * Math.PI, 0.5 * D)) / ((1 - cpd.opt.outliers) * N);
        if (outlier_tmp == 0) outlier_tmp = 1e-10;

        for (n = 0; n < N; n++) {
            sp = 0;
            for (m = 0; m < M; m++) {
                razn = 0;
                for (d = 0; d < D; d++) {
                    diff = cpd.X[n][d] - cpd.registeredY[m][d];
                    diff *= diff;
                    razn += diff;
                }
                double value = Math.exp(razn / ksig);
                P[m] = value;
                sp += value;
            }
            sp += outlier_tmp;
            for (m = 0; m < M; m++) {
                temp_x = P[m] / sp;

                if (n == 0) {
                    P1[m] = temp_x;
                    C[m] = n + 1;
                }
                ;

                if (temp_x > P1[m]) {
                    P1[m] = temp_x;
                    C[m] = n + 1;
                }
            }
        }
        return C;
    }


    /**
     * The non-rigid CPD point-set registration
     *
     * @param cpd {@link CPD}
     * @return registered Y point set
     */
    private static double[][] cpd_GRBF(CPD cpd) {
        double lambda = cpd.opt.lambda;
        double beta = cpd.opt.beta;
        double tol = cpd.opt.tol;
        double outlier = cpd.opt.outliers;
        int max_it = cpd.opt.max_iter;

        //int N = rowCount(cpd.X);
        int D = colCount(cpd.X);
        int M = rowCount(cpd.Y);

        //Initialization
        double[][] T = copy(cpd.Y);
        int iter = 0;
        double ntol = tol + 10;
        double[][] W = new double[M][D];

        double sigma2 = getSigma2(cpd);
        double sigma2_init = sigma2;
        //Construct affinity matrix G
        double[][] G = cpd_G(cpd.Y, cpd.Y, beta);
        double L = 1;
        double L_old;

        Plot plot = null;
        if (cpd.opt.plot) {
            plot = new Plot();
            plot.show();
            plot.addXYSeries("X", cpd.X);
            plot.addXYSeries("T", T);
        }

        while ((iter < max_it) && (ntol > tol) && (sigma2 > 1e-8)) //(sigma2 > 1e-8)
        {
            //I don't implement the Fast Gauss Transform (FGT)
            L_old = L;
            ProbabilityMatrix P;
            if (cpd.opt.fgt) {
                sigma2 = Math.max(0.05, sigma2);
                P = Utils.cpd_P_FGT(cpd.X, T, sigma2, outlier, sigma2_init);
                if (cpd.opt.print)
                    System.out.print("(fgt) ");
            } else
                P = Utils.cpd_P(cpd.X, T, sigma2, outlier);
            L = P.E + lambda / 2.0 * ArraysMath.trace(ArraysMath.mmul(ArraysMath.mmul(transpose(W), G), W));
            ntol = Math.abs((L - L_old) / L);

            if (cpd.opt.print)
                System.out.println(
                        "CPD non-rigid" + " : dL= " + ntol +
                                ", iter= " + iter + ", sigma2= " + sigma2);

            double[][] W_temp1 = new double[M][M];
            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                    W_temp1[i][j] = G[i][j] * P.P1[i];
            for (int i = 0; i < M; i++)
                W_temp1[i][i] += lambda * sigma2;

            double[][] W_temp2 = new double[M][D];
            for (int i = 0; i < M; i++)
                for (int j = 0; j < D; j++)
                    W_temp2[i][j] = P.PX[i][j] - cpd.Y[i][j] * P.P1[i];


            W = ArraysMath.solve(W_temp1, W_temp2);
            double[][] GW = ArraysMath.mmul(G, W);
            for (int i = 0; i < T.length; i++)
                for (int j = 0; j < T[0].length; j++)
                    T[i][j] = cpd.Y[i][j] + GW[i][j];

            if (cpd.opt.plot) {
                plot.setXYSeries("T", T);
                try {
                    Thread.sleep(60);
                } catch (InterruptedException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }

            double Np = ArraysMath.sum(P.P1);
            double temp1 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(cpd.X, cpd.X), P.Pt1, ArraysMath.col)));
            double temp2 = ArraysMath.sum(ArraysMath.sum(ArraysMath.mul(ArraysMath.mul(T, T), P.P1, ArraysMath.col)));
            double temp3 = -2 * ArraysMath.trace(ArraysMath.mmul(transpose(P.PX), T));
            cpd.sigma2 = sigma2;
            sigma2 = Math.abs((temp1 + temp2 + temp3) / (Np * D));
            iter++;
        }
        System.out.println("CPD registration succesfully completed.");
        return T;
    }

    protected static double getSigma2(CPD cpd) {
        int N = rowCount(cpd.X);
        int D = colCount(cpd.X);
        int M = rowCount(cpd.Y);
        double sigma2 =
                (M * ArraysMath.trace(ArraysMath.mmul(cpd.X, transpose(cpd.X))) +
                        N * ArraysMath.trace(ArraysMath.mmul(cpd.Y, transpose(cpd.Y))) -
                        2.0 * ArraysMath.dot(ArraysMath.sum(cpd.X), ArraysMath.sum(cpd.Y)))
                        / (M * D * M);
        return sigma2;
    }

}
