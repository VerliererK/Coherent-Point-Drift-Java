package com.cklin.cpd;

import com.cklin.math.ArraysMath;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import static com.cklin.math.ArraysUtils.*;

public class NonRigid {
    /**
     * The non-rigid CPD point-set registration
     *
     * @param cpd {@link CPD}
     * @return registered Y point set
     */
    public static double[][] GRBF(CPD cpd) {
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

        double sigma2 = Utils.getSigma2(cpd);
        double sigma2_init = sigma2;
        //Construct affinity matrix G
        double[][] G = Utils.cpd_G(cpd.Y, cpd.Y, beta);
        double L = 1;
        double L_old;

        while ((iter < max_it) && (ntol > tol) && (sigma2 > 1e-8)) //(sigma2 > 1e-8)
        {
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

    public static double[][] GRBF_lowrank(CPD cpd) {
        double lambda = cpd.opt.lambda;
        double beta = cpd.opt.beta;
        double tol = cpd.opt.tol;
        double outlier = cpd.opt.outliers;
        int numeig = cpd.opt.numeig;
        int max_it = cpd.opt.max_iter;

        //int N = rowCount(cpd.X);
        int D = colCount(cpd.X);
        int M = rowCount(cpd.Y);
        if (numeig < 1)
            numeig = (int) Math.sqrt(M);
        //Initialization
        double[][] T = copy(cpd.Y);
        int iter = 0;
        double ntol = tol + 10;
        double L = 1;
        double L_old;
        double[][] W = new double[M][D];

        double sigma2 = Utils.getSigma2(cpd);
        double sigma2_init = sigma2;
        //Find 'numeig' eigenvectors of the affinity matrix G
        //TODO:
        double[][] G = Utils.cpd_G(cpd.Y, cpd.Y, beta);
        //[Q,S]=eigs(G,numeig,'lm',OPTS);
        //EigenDecomposition eig = new EigenDecomposition(G);
        double[][] Q = null;
        double[][] S = null;
        //[Q,S] = GRBF_lowrankQS(cpd.Y, beta, numeig);
        //invS=spdiags(1./diag(S),0,numeig,numeig);

        while ((iter < max_it) && (ntol > tol) && (sigma2 > 1e-8)) //(sigma2 > 1e-8)
        {
            L_old = L;
            double[][] QtW = ArraysMath.mmul(transpose(Q), W);

            ProbabilityMatrix P;
            if (cpd.opt.fgt) {
                sigma2 = Math.max(0.05, sigma2);
                P = Utils.cpd_P_FGT(cpd.X, T, sigma2, outlier, sigma2_init);
                if (cpd.opt.print)
                    System.out.print("(fgt) ");
            } else
                P = Utils.cpd_P(cpd.X, T, sigma2, outlier);

            L = P.E + lambda / 2.0 * ArraysMath.trace(ArraysMath.mmul(ArraysMath.mmul(transpose(QtW), S), QtW));
            ntol = Math.abs((L - L_old) / L);

            if (cpd.opt.print)
                System.out.println(
                        "CPD non-rigid (lowrank)" + " : dL= " + ntol +
                                ", iter= " + iter + ", sigma2= " + sigma2);

            //M-step. Solve linear system for W.
            //TODO:
            // dPQ=dP*Q;
            // F=PX-dP*Y;
            // W=1/(lambda*sigma2)*(F-dPQ*  ( (lambda*sigma2*invS+Q'*dPQ) \ (Q'*F) ) );

            // update Y postions
            // same as T=Y+G*W for full rank;
            //T=Y+(Q*(S*(Q'*W)));
            double[][] QSQW = ArraysMath.mmul(Q, ArraysMath.mmul(S, ArraysMath.mmul(transpose(Q), W)));
            for (int i = 0; i < T.length; i++)
                for (int j = 0; j < T[0].length; j++)
                    T[i][j] = cpd.Y[i][j] + QSQW[i][j];

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

    private static void GRBF_lowrankQS() {
        //TODO:
    }
}
