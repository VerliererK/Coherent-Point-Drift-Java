package test;

import java.awt.Color;
import java.awt.Shape;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.UnknownKeyException;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


public class Plot {
    private XYSeriesCollection dataset = new XYSeriesCollection();

    public XYSeries getseries(String name) {
        try {
            return dataset.getSeries(name);
        } catch (UnknownKeyException e) {
            return null;
        }
    }

    public void addXYSeries(String name, double[][] data) {
        if (getseries(name) != null)
            return;
        XYSeries s = new XYSeries(name, false);
        for (int i = 0; i < data.length; i++) {
            double x = data[i][0];
            double y = data[i][1];
            s.add(x, y);
        }
        dataset.addSeries(s);
    }

	/*public void addXYSeries(String name, RealMatrix data){
        if(getseries(name) != null)
			return;
		XYSeries s = new XYSeries(name, false);
		for (int i = 0; i < data.getRowDimension(); i++) {
			double x = data.getEntry(i, 0);
			double y = data.getEntry(i, 1);
			s.add(x, y);
		}
		dataset.addSeries(s);
	}*/

    public void setXYSeries(String name, double[][] data) {
        XYSeries s = getseries(name);
        if (s == null)
            return;
        s.clear();
        for (int i = 0; i < data.length; i++) {
            double x = data[i][0];
            double y = data[i][1];
            s.add(x, y);
        }
    }
    /*public void setXYSeries(String name, RealMatrix data){
        XYSeries s = getseries(name);
		if(s == null)
			return;
		s.clear();
		for (int i = 0; i < data.getRowDimension(); i++) {
			double x = data.getEntry(i, 0);
			double y = data.getEntry(i, 1);
			s.add(x, y);
		}
	}*/

    private JFreeChart chart = ChartFactory.createXYLineChart("", "", "", dataset);

    public void text(String title, String xAxisLabel, String yAxisLabel) {
        chart.setTitle(title);
        chart.getXYPlot().getDomainAxis().setLabel(xAxisLabel);
        chart.getXYPlot().getRangeAxis().setLabel(yAxisLabel);
    }

    private XYLineAndShapeRenderer render = new XYLineAndShapeRenderer(false, true);

    public void render(String series, Color color, Shape shape) {
        int s = dataset.getSeriesIndex(series);
        render.setSeriesPaint(s, color);
        render.setSeriesShape(s, shape);
        ((XYPlot) chart.getPlot()).setRenderer(render);
    }

    private void ini_render() {
        chart.getPlot().setBackgroundPaint(Color.white);
        render.setBaseShapesFilled(false);
    }

    private jf figure = new jf("figure", "figure", chart, 500, 500);

    public void show() {
        ini_render();
        ((XYPlot) chart.getPlot()).setRenderer(render);
        figure.pack();
        figure.setVisible(true);
    }

    public void figure_title(String title) {
        figure.setTitle(title);
    }

    public void figure_size(int width, int height) {
        figure.setPreferredSize(new java.awt.Dimension(width, height));
    }

    private class jf extends JFrame {
        /**
         *
         */
        private static final long serialVersionUID = 1L;

        public jf(String applicationTitle, String chartTitle, JFreeChart chart, int width, int height) {
            super(applicationTitle);
            // This will create the dataset
            //PieDataset dataset = createDataset();
            // based on the dataset we create the chart
            //JFreeChart chart = createChart(dataset, chartTitle);
            // we put the chart into a panel
            ChartPanel chartPanel = new ChartPanel(chart);
            // default size
            chartPanel.setPreferredSize(new java.awt.Dimension(width, height));
            // add it to our application
            setContentPane(chartPanel);
        }
    }
}
