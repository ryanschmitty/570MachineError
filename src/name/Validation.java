package name;

import java.util.*;
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

public class Validation {

	/**
	 * @param args
	 * @throws DriverException
	 * @throws SQLException
	 */
	public static void main(String[] args) throws Exception {
		CPLOPConnection conn = new CPLOPConnection();

		FileWriter outputFile = new FileWriter("output.csv");
		BufferedWriter output = new BufferedWriter(outputFile);
		testCorrelation(conn, output);
		testDiscorrelation(conn, output);
		output.close();
		outputFile.close();

		/*
		ArrayList<Double> l1 = new ArrayList<Double>(); 
		l1.add(0,1.0);
		l1.add(1,2.0); 
		l1.add(2,2.0); 
		l1.add(3,1.0); 
		l1.add(4,3.0);
		 
		
		ArrayList<Double> l2 = new ArrayList<Double>(); 
		l2.add(0,3.0);
		l2.add(1,4.0); 
		l2.add(2,2.0); 
		l2.add(3,1.0); 
		l2.add(4,2.0); 
		
		double d = mahalanobisCorrelation(l1, l2, 5);
		System.out.println("MahalanobisDistance = " + d);*/
		
	}

	private static void testCorrelation(CPLOPConnection conn,
			BufferedWriter output) throws SQLException, IOException {
		List<Map<String, Object>> isolates = conn.getIsolates();
		List<List<Double>> histograms = new ArrayList<List<Double>>();

		int count = 0;
		for (Map<String, Object> isolateMetaDatas : isolates) {
			// System.out.println(isolateMetaDatas);
			if ((Integer) isolateMetaDatas.get("count") < 5) {
				break;
			}
			count++;
			List<Map<String, Object>> pyroprints = conn
					.getPyroprints((String) isolateMetaDatas.get("isolate"));

			int pyroprintCount = 0;
			Map<String, Object> pyroprintMetadataModel = pyroprints.get(0);
			for (Map<String, Object> pyroprintMetadata : pyroprints) {
				// System.out.println(pyroprintMetadata);

				// pyroprintMetadata = pyroprintMetadataModel;
				if (equals(pyroprintMetadataModel, pyroprintMetadata)) {
					pyroprintCount++;
					Integer pyroId = (Integer) pyroprintMetadata
							.get("pyroprint");
					List<Map<String, Object>> histogram = conn
							.getHistogram(pyroId);
					// System.out.println(histogram);
					ArrayList<Double> peakHeights = new ArrayList<Double>();
					for (Map<String, Object> tuple : histogram) {
						peakHeights.add((Double) tuple.get("peakHeight"));
					}

					histograms.add(peakHeights);

				}
			}

			validate(output, (String) isolateMetaDatas.get("isolate") + " "
					+ pyroprintMetadataModel.get("region") + ","
					+ pyroprintCount, histograms);
		}
	}

	private static void testDiscorrelation(CPLOPConnection conn,
			BufferedWriter output) throws SQLException, IOException {
		List<Map<String, Object>> isolates = conn.getIsolates();
		List<List<Double>> histograms = new ArrayList<List<Double>>();

		int count = 0;
		for (Map<String, Object> isolateMetaDatas : isolates) {
			// System.out.println(isolateMetaDatas);
			if ((Integer) isolateMetaDatas.get("count") < 5) {
				break;
			}
			count++;
			List<Map<String, Object>> pyroprints = conn
					.getPyroprints((String) isolateMetaDatas.get("isolate"));

			Map<String, Object> pyroprintMetadata = pyroprints.get(0);

			Integer pyroId = (Integer) pyroprintMetadata.get("pyroprint");
			List<Map<String, Object>> histogram = conn.getHistogram(pyroId);
			// System.out.println(histogram);
			ArrayList<Double> peakHeights = new ArrayList<Double>();
			for (Map<String, Object> tuple : histogram) {
				peakHeights.add((Double) tuple.get("peakHeight"));
			}

			histograms.add(peakHeights);
		}

		validate(output, "Discorrelation," + count + ",", histograms);

	}

	private static boolean equals(Map<String, Object> model,
			Map<String, Object> other) {
		return model.get("sequencePrimer").equals(other.get("sequencePrimer"))
				&& model.get("region").equals(other.get("region"))
				&& model.get("forwardPrimer")
						.equals(other.get("forwardPrimer"))
				&& model.get("forwardPrimer")
						.equals(other.get("forwardPrimer"))
				&& model.get("dispensation").equals(other.get("dispensation"));
	}

	private static void validate(BufferedWriter output, String name,
			List<List<Double>> histograms) throws IOException {
		int size = histograms.size();
		Double[][] pearsonTable = new Double[size][size];
		Double[][] jaccardTable = new Double[size][size];
		Double[][] mahalanobisTable = new Double[size][size];
		Double mAverage = 0.0, jAverage = 0.0, pAverage = 0.0;

		for (int range = 60; range <= 104; range++) {
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					pearsonTable[i][j] = pearsonCorrelation(histograms.get(i),
							histograms.get(j), range);
				}
			}
			pAverage = getAverage(pearsonTable);

			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					jaccardTable[i][j] = jaccardCorrelation(histograms.get(i),
							histograms.get(j), range);
				}
			}
			jAverage = getAverage(jaccardTable);

			
			  for (int i = 0; i < size; i++) { 
				  for (int j = 0; j < size; j++) {
					  mahalanobisTable[i][j] = mahalanobisCorrelation(histograms
							  .get(i), histograms.get(j), range); 
				  } 
			  } 
			  mAverage = getAverage(mahalanobisTable);
			 

			output.write(name + "," + range + "," + pAverage + "," + jAverage
					+ "," + mAverage + "\n");
		}
	}

	private static Double getAverage(Double[][] table) {
		Double sum = 0.0;
		for (int i = 0; i < table.length; i++) {
			for (int j = 0; j < table.length; j++) {
				sum += table[i][j];
			}
		}
		return sum / (table.length * table.length);
	}

	/**
	 * http://en.wikipedia.org/wiki/Pearson_product-
	 * moment_correlation_coefficient
	 * 
	 * @param firstPyroprint
	 * @param secondPyroprint
	 * @param range
	 * @return
	 */
	public static Double pearsonCorrelation(List<Double> firstPyroprint,
			List<Double> secondPyroprint, int range) {
		Double firstAverage = 0.0;
		Double secondAverage = 0.0;
		Double firstDiff = 0.0;
		Double secondDiff = 0.0;
		Double topSum = 0.0;

		firstAverage = listAverage(firstPyroprint, range);
		secondAverage = listAverage(secondPyroprint, range);

		for (int i = 0; i < Math.min(firstPyroprint.size(), secondPyroprint
				.size())
				&& i < range; i++) {
			firstDiff += Math.pow(firstPyroprint.get(i) - firstAverage, 2);
			secondDiff += Math.pow(secondPyroprint.get(i) - secondAverage, 2);
			topSum += (firstPyroprint.get(i) - firstAverage)
					* (secondPyroprint.get(i) - secondAverage);
		}

		return (Double) topSum / (Math.sqrt(firstDiff) * Math.sqrt(secondDiff));

	}

	private static Double listAverage(List<Double> list, int range) {
		Double sum = 0.0;
		for (int i = 0; i < list.size() && i < range; i++) {
			sum += list.get(i);
		}
		return sum / list.size();
	}

	private static Double mahalanobisCorrelation(List<Double> firstPyroprint,
			List<Double> secondPyroprint, int range) {
		double[] first = listToArray(firstPyroprint, range);
		double[] sec = listToArray(secondPyroprint, range);
		double firstMean = 0.0;
		double secMean = 0.0;
		for (int i = 0; i < range; i++) {
			firstMean += first[i];
			secMean = sec[i];
		}
		firstMean /= range;
		secMean /= range;

		MahalanobisDistanceMeasure mah = new MahalanobisDistanceMeasure();
		DenseVector v1 = new DenseVector(first);
		DenseVector v2 = new DenseVector(sec);
		double[] meanV = new double[2];
		meanV[0] = firstMean;
		meanV[1] = secMean;
		
		double[][] arr = new double[range][2];
		for(int x = 0; x < range; x++)
		{
			arr[x][0]= first[x];
			arr[x][1]= sec[x];
		}
		/*double[][] arr = new double[range][2];
        arr[0] = first;
        arr[1] = sec;*/
		for(int i = 0; i < range; i++)
		{
			System.out.println(arr[i][0] + ", " + arr[i][1]);
		}
		
		
		
		Covariance cov = new Covariance(arr);
		RealMatrix m = cov.getCovarianceMatrix();
		System.out.println("--- " + m.toString());
/*		double[][] m = new double[5][5];
		m[0][0] = 1;
		m[1][1] = 1;
		m[2][2] = 1;
		m[3][3] = 1;
		m[4][4] = 1;*/
		
		DenseMatrix covMatrix = new DenseMatrix(m.getData());
		DenseVector meanVector = new DenseVector(meanV);
		mah.setMeanVector(meanVector);
		mah.setCovarianceMatrix(covMatrix);
		return mah.distance(v1, v2);
	}

	private static Double jaccardCorrelation(List<Double> firstPyroprint,
			List<Double> secondPyroprint, int range) {
		double[] first = listToArray(firstPyroprint, range);
		double[] sec = listToArray(secondPyroprint, range);

		TanimotoDistanceMeasure dist = new TanimotoDistanceMeasure();
		DenseVector v1 = new DenseVector(first);
		DenseVector v2 = new DenseVector(sec);

		return 2 * (1 - dist.distance(v1, v2)) - 1;
	}

	private static double[] listToArray(List<Double> list, int range) {
		double[] arr = new double[range];
		for (int i = 0; i < range && i < list.size(); i++) {
			arr[i] = list.get(i);
		}
		return arr;
	}
}
