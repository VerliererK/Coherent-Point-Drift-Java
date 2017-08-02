package com.cklin.cpd;

/**
 * Structure that can be used to re-scale and shift the point
 */
class Normal {
    /**
     * shift variables
     */
    protected double[] xd, yd;
    /**
     * re-scale variables
     */
    protected double xScale, yScale;

    /**
     * Structure that can be used to re-scale and shift the point
     *
     * @param xd     shift variables
     * @param yd     shift variables
     * @param xScale re-scale variables
     * @param yScale re-scale variables
     */
    protected Normal(double[] xd, double[] yd, double xScale, double yScale) {
        this.xd = xd.clone();
        this.yd = yd.clone();
        this.xScale = xScale;
        this.yScale = yScale;
    }
}
