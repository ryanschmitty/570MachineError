package name;

import java.util.*;
import java.lang.Math;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import org.apache.mahout.common.distance.MahalanobisDistanceMeasure;
import org.apache.mahout.common.distance.TanimotoDistanceMeasure;
import org.apache.mahout.math.*;
import org.apache.mahout.math.Vector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.correlation.Covariance;

import name.CPLOPConnection.DriverException;

public class DetectMachineError {

    /**
     * @param args
     * @throws DriverException
     * @throws SQLException
     */
    public static void main(String[] args) throws Exception {
        CPLOPConnection conn = new CPLOPConnection();

//        List<Double> comparisonVals = getPearsonComparisonValues(2377, 2377, conn);
//
//        System.out.println("" + comparisonVals.size() + " Comparison Values:");
//        Iterator it = comparisonVals.iterator();
//        while (it.hasNext()) {
//            System.out.println(it.next());
//        }
//        System.out.println("Done!");
        
        //Look for a 15% change in correlation between two lengths of comparison
//        Double threshold = 0.15;
//        Integer len = hasCorrelationDrop(2377, 2377, threshold, conn);
//        if (len < 0)
//            System.out.println("No correlation change greater than threshold.");
//        else
//            System.out.println("Correlation changed more than "+threshold+" between lengths "+(len-1)+" and "+len+".");

        //Test average and standard deviation calculation
//        ArrayList<Double> testList = new ArrayList<Double>();
//        testList.add(2.0);
//        testList.add(4.0);
//        testList.add(4.0);
//        testList.add(4.0);
//        testList.add(5.0);
//        testList.add(5.0);
//        testList.add(7.0);
//        testList.add(9.0);
//        System.out.println("List Values: "+testList);
//        System.out.println("Average: "+calculateAverage(testList));
//        System.out.println("Std Dev: "+calculateStdDev(testList, calculateAverage(testList)));

        //Check for deviations from the standard deviation
//        Double numStandardDeviations = 3.0;
//        Integer windowSize = 10;
//        List<Double> inputPeakHeights = getPeakHeightHistogram(2377, conn);
//        List<Double> comparisonPeakHeights = getPeakHeightHistogram(2377, conn);
//        Integer index = getPeakHeightOutlierIndex(inputPeakHeights, comparisonPeakHeights, numStandardDeviations, windowSize);
//        if (index < 0)
//            System.out.println("No outlier peak height detected.");
//        else
//            System.out.println("Outlier peak height detected.");

        //Check for Pearson correlation change greater than some number of
        //running standard deviations from the running average
        Double numStandardDeviations = 2.0;
        Integer windowSize = 10;
        Integer index = hasCorrelationDropWithRunningStdDev(2377, 2377, numStandardDeviations, windowSize, conn);
        if (index < 0)
            System.out.println("No Pearson Correlation anomaly detected.");
        else
            System.out.println("Outlier Pearson Correlation drop detected at length "+index+".");
    }

    /**
     * Check for a drop in correlation and report length at which it occured.
     * @param pyroId1 The first pyroprint
     * @param pyroId2 The second pyroprint
     * @param correlationChangeThreshold The upper-bound on accetable correlation change
     * @param conn The CPLOPConnection
     * @return The length value of the first correlation change that is greater than the
     * provided threshold, or -1 if no correlation change is greater than the threshold.
     */
    public static Integer hasCorrelationDrop(Integer pyroId1, Integer pyroId2, Double correlationChangeThreshold, CPLOPConnection conn) throws SQLException {
        //Get a list of correlations for every possible length (e.g. 1..104)
        List<Double> comparisonVals = getPearsonComparisonValues(pyroId1, pyroId2, conn);

        //Compare all correlation values with the previous one and report if the
        //change is greater than the provided threshold.
        Iterator<Double> it = comparisonVals.iterator();
        Double prev = 0.0;
        Integer cur = 0;
        if (it.hasNext()) prev = it.next();
        while (it.hasNext()) {
            cur++;
            // No absolute value so that an increase in correlation doesn't
            // count.
            if (prev - it.next() > correlationChangeThreshold) return cur;
        }
        return -1; //No change was greater than the threshold
    }

    /**
     * Search for the first index where the Pearson Correlation falls more than
     * a specified number of standard deviations.
     * @param pyroId1 The first pyroprint.
     * @param pyroId2 The second pyroprint.
     * @param numStdDevs The number of acceptable standard deviations.
     * @param windowSize The window size to use in running average and standard
     * deviation computation.
     * @param conn The CPLOPConnection
     * @return The length of the first correlation drop that was greater than
     * the specified number of standard deviations, or -1 if no correlation
     * change was outisde the acceptable standard deviations.
     */
    public static Integer hasCorrelationDropWithRunningStdDev(Integer pyroId1, Integer pyroId2, Double numStdDevs, Integer windowSize, CPLOPConnection conn) throws SQLException {
        //Get a list of correlations for every possible length (e.g. 1..104)
        List<Double> pearsonVals = getPearsonComparisonValues(pyroId1, pyroId2, conn);

        ListIterator<Double> inputIterator = pearsonVals.listIterator(windowSize); //start at index=windowSize
        for (int i=windowSize; i<pearsonVals.size(); ++i) {
            List<Double> window = pearsonVals.subList(i-windowSize, i);
            Double avg = calculateAverage(window); //average Pearson Correlation over the window
            Double stdDev = calculateStdDev(window, avg); //standard devation of Pearson Correlation over the window
            //Return index if its Pearson value falls more than the specified
            //number of standard deviations.
            if (avg - pearsonVals.get(i) > numStdDevs*stdDev) {
                System.out.println("Pearson Correlation value ("+pearsonVals.get(i)+") greater than "+numStdDevs+" standard deviations ("+stdDev+") from running average ("+avg+") of length "+windowSize+" at index "+i+".");
                return i;
            }

        }
        return -1; //No correlation drop large enough
    }

    /**
     * Search for first peak height outlier, in terms of standard deviation.
     * @param inputPeakHeights The histogram of peak heights you want to examine
     * @param comparisonPeakHeights The histogram of good or average peak
     * heights you want to compare against.
     * @param numStdDevs The number of acceptable standard deviations
     * @param windowSize The window size to use in running average and standard
     * deviation computation.
     * @return The index in inputPeakHeights where its running average is
     * outside of the provided standard devation range, or -1 if it never is
     * ouside the range.
     */
    public static Integer getPeakHeightOutlierIndex(List<Double> inputPeakHeights, List<Double> comparisonPeakHeights, Double numStdDevs, Integer windowSize) {
        //Error if histograms aren't long enough, or erroneous windowSize
        if (inputPeakHeights.size() <= windowSize || comparisonPeakHeights.size() <= windowSize || windowSize < 1) return -1;

        //Iterate over inputPeakHeights from [windowSize, size) and check for
        //deviation from the running average of the comparison values by the
        //specified number of standard deviations.
        ListIterator<Double> inputIterator = inputPeakHeights.listIterator(windowSize);
        Integer len = inputPeakHeights.size() < comparisonPeakHeights.size() ? inputPeakHeights.size() : comparisonPeakHeights.size();
        for (int i=windowSize; i<len; ++i) {
            //Get the running average from the input values.
            Double inputAvg = calculateAverage(comparisonPeakHeights.subList(i-windowSize, i));
            //Get the running average and standard devation from
            //the comparison values.
            List<Double> comparisonSubList = comparisonPeakHeights.subList(i-windowSize, i);
            Double comparisonAvg = calculateAverage(comparisonSubList);
            Double comparisonStdDev = calculateStdDev(comparisonSubList, comparisonAvg);
            //Check for devation from the average by the specified number of
            //standard devations.
            if (Math.abs(comparisonAvg - inputAvg) > numStdDevs*comparisonStdDev) {
                System.out.println("Peak height average ("+inputAvg+") greater than "+numStdDevs+" standard deviations ("+comparisonStdDev+") from running average ("+comparisonAvg+") of length "+windowSize+" at index "+i+".");
                return i;
            }
        }

        return -1; //No value deviated enough to be flagged
    }

    /**
     * Compares two pyroprints using the Pearson Correlation at varying lengths.
     * @return A List of comparison values for each length from 1..length.
     */
    public static List<Double> getPearsonComparisonValues(Integer pyroId1, Integer pyroId2, CPLOPConnection conn) throws SQLException {
        //Variable init
        ArrayList<Double> comparisonValues = new ArrayList<Double>();

        //Get Peak Height Histograms from Pyro IDs
        List<Double> peakHeights1 = getPeakHeightHistogram(pyroId1, conn);
        List<Double> peakHeights2 = getPeakHeightHistogram(pyroId2, conn);

        //Only compare up to shortest length
        int len = peakHeights1.size() < peakHeights2.size() ? peakHeights1.size() : peakHeights2.size();

        //Add comparison values for all possible lengths
        for (int i=1; i<len; ++i) {
            comparisonValues.add( Validation.pearsonCorrelation(peakHeights1, peakHeights2, i) );
        }

        return comparisonValues;
    }

    /**
     * Simple utility function to get the Peak Height Histogram from a Pyro ID.
     */
    private static List<Double> getPeakHeightHistogram(Integer pyroId, CPLOPConnection conn) throws SQLException {
        List<Map<String, Object>> histogram = conn.getHistogram(pyroId);
        ArrayList<Double> peakHeights = new ArrayList<Double>();
        for (Map<String, Object> tuple : histogram) {
            peakHeights.add((Double)tuple.get("peakHeight"));
        }
        return peakHeights;
    }

    /**
     * Simple average computation helper.
     */
    private static Double calculateAverage(List<Double> values) {
        Double average = 0.0;
        Iterator<Double> it = values.iterator();
        while (it.hasNext()) {
            average += it.next();
        }
        return average / values.size();
    }

    /**
     * Simple utility funciton to compute the standard deviation of a set of
     * values.
     * @param values The list of double precision values
     * @param average The precomputed average of the values (so we don't do lots
     * of extra work).
     * @return The standard devation of the list of values; note that this will
     * be incorrect if you do not provide the correct average.
     */
    private static Double calculateStdDev(List<Double> values, Double average) {
        Double stdDeviation = 0.0;
        Iterator<Double> it = values.iterator();
        while (it.hasNext()) {
            Double oneValue = it.next() - average;
            stdDeviation += oneValue * oneValue;
        }
        stdDeviation = Math.sqrt(stdDeviation/values.size());
        return stdDeviation;
    }
}
