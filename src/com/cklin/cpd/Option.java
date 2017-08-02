package com.cklin.cpd;

/**
 * Set of CPD options
 */
public class Option {
    /**
     * Only implement Non-rigid
     */
    public Method method = Method.NonRigid;
    /**
     * the width of Gaussian kernel (smoothness)<br>
     * strictly non-rigid params
     */
    public double beta = 2.0;
    /**
     * regularization weight<br>
     * strictly non-rigid params
     */
    public double lambda = 3.0;
    /**
     * noise weight
     */
    public double outliers = 0.4;
    /**
     * max number of iterations
     */
    public int max_iter = 100;
    /**
     * tolerance
     */
    public double tol = 1e-10;
    /**
     * The number of largest eigenvectors to use,
     * </br>try sqrt(N) sqrt(M), where M is the length of Y.
     */
    public int numeig = 30;
    /**
     * Use a Fast Gauss transform (FGT)
     */
    public boolean fgt = false;
    /**
     * normalize to unit variance and zero mean before registering (default = true)
     */
    public boolean normalize = true;
    /**
     * compute correspondence vector at the end of registration (not being estimated by default)
     */
    public boolean corresp = false;
    /**
     * show every iteration
     */
    public boolean print = true;
    /**
     * show every iteration with {@link Plot}
     */
    public boolean plot = false;
    /**
     * For rigid: true- strict rotation, false- allow reflections.
     */
    public boolean rot = true;
    /**
     * For rigid: true- estimate scaling, false- don't estimate scaling.
     */
    public boolean scale = true;
}
