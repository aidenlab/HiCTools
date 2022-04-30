/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */


package hic.tools.utils.original;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.expected.ExpectedValueFunctionImpl;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Computes an "expected" density vector.  Essentially there are 3 steps to using this class
 * <p/>
 * (1) instantiate it with a collection of Chromosomes (representing a genome) and a grid size
 * (2) loop through the pair data,  calling addDistance for each pair, to accumulate all counts
 * (3) when data loop is complete, call computeDensity to do the calculation
 * <p/>
 * <p/>
 * Methods are provided to save the result of the calculation to a binary file, and restore it.  See the
 * DensityUtil class for example usage.
 *
 * @author Jim Robinson
 * @since 11/27/11
 */
public class ExpectedValueCalculation {

    private final int gridSize;

    private final int numberOfBins;
    /**
     * Map of chromosome index -> total count for that chromosome
     */
    private final Map<Integer, Double> chromosomeCounts = new ConcurrentHashMap<>();
    /**
     * Map of chromosome index -> "normalization factor", essentially a fudge factor to make
     * the "expected total"  == observed total
     */
    private final Map<Integer, Double> chrScaleFactors = new ConcurrentHashMap<>();
    private final NormalizationType type;
    // A little redundant, for clarity
	/**
	 * Genome wide count of binned reads at a given distance
	 */
	private final double[] actualDistances;
	/**
	 * Expected count at a given binned distance from diagonal
	 */
	private ListOfDoubleArrays densityAvg;
	/**
	 * Chromosome in this genome, needed for normalizations
	 */
	private final Map<Integer, Chromosome> chromosomesMap = new ConcurrentHashMap<>();

    /**
     * Instantiate a DensityCalculation.  This constructor is used to compute the "expected" density from pair data.
     *
     * @param chromosomeHandler Handler for list of chromosomesMap, mainly used for size
     * @param gridSize          Grid size, used for binning appropriately
     * @param type              Identifies the observed matrix type,  either NONE (observed), VC, or KR.
     */
    public ExpectedValueCalculation(ChromosomeHandler chromosomeHandler, int gridSize, NormalizationType type) {

        this.type = type;
        this.gridSize = gridSize;

        long maxLen = 0;

        for (Chromosome chr : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {
            if (chr != null) {
                chromosomesMap.put(chr.getIndex(), chr);
                try {
                    maxLen = Math.max(maxLen, chr.getLength());
                }
                catch (NullPointerException error) {
                    System.err.println("Problem with creating fragment-delimited maps, NullPointerException.\n" +
                            "This could be due to a null fragment map or to a mismatch in the chromosome name in " +
                            "the fragment map vis-a-vis the input file or chrom.sizes file.\n" +
                            "Exiting.");
                    System.exit(63);
                }
                catch (ArrayIndexOutOfBoundsException error) {
                    System.err.println("Problem with creating fragment-delimited maps, ArrayIndexOutOfBoundsException.\n" +
                            "This could be due to a null fragment map or to a mismatch in the chromosome name in " +
                            "the fragment map vis-a-vis the input file or chrom.sizes file.\n" +
                            "Exiting.");
                    System.exit(22);
                }
            }
        }

        numberOfBins = (int) (maxLen / gridSize) + 1;

        actualDistances = new double[numberOfBins];
        Arrays.fill(actualDistances, 0);
    }

    public int getGridSize() {
        return gridSize;
    }


    /**
     * Add an observed distance.  This is called for each pair in the data set
     *
     * @param chrIdx index of chromosome where observed, so can increment count
     * @param bin1   Position1 observed in units of "bins"
     * @param bin2   Position2 observed in units of "bins"
     */
    public synchronized void addDistance(Integer chrIdx, int bin1, int bin2, double weight) {

        // Ignore NaN values    TODO -- is this the right thing to do?
        if (Double.isNaN(weight)) return;

        int dist;
        Chromosome chr = chromosomesMap.get(chrIdx);
        if (chr == null) return;

        if (chromosomeCounts.containsKey(chrIdx)) {
            double count = chromosomeCounts.get(chrIdx);
            chromosomeCounts.put(chrIdx, count + weight);
        } else {
            chromosomeCounts.put(chrIdx, weight);
        }
        dist = Math.abs(bin1 - bin2);

        actualDistances[dist] += weight; // Math.log(1 + weight);
    }

    public void merge(ExpectedValueCalculation otherEVCalc) {
        for (Map.Entry<Integer, Chromosome> entry : otherEVCalc.chromosomesMap.entrySet()) {
            Chromosome chr = chromosomesMap.get(entry.getKey());
            if (chr != null) {
                if (otherEVCalc.chromosomeCounts.get(entry.getKey()) != null) {
                    Double count = chromosomeCounts.get(entry.getKey());
                    if (count == null) {
                        chromosomeCounts.put(entry.getKey(), otherEVCalc.chromosomeCounts.get(entry.getKey()));
                    } else {
                        chromosomeCounts.put(entry.getKey(), count + otherEVCalc.chromosomeCounts.get(entry.getKey()));
                    }
                }
            }
        }
        for (int i = 0; i < actualDistances.length; i++) {
            actualDistances[i] += otherEVCalc.actualDistances[i];
        }
    }

    public boolean hasData() {
        return !chromosomeCounts.isEmpty();
    }

    /**
     * Compute the "density" -- port of python function getDensityControls().
     * The density is a measure of the average distribution of counts genome-wide for a ligated molecule.
     * The density will decrease as distance from the center diagonal increases.
     * First compute "possible distances" for each bin.
     * "possible distances" provides a way to normalize the counts. Basically it's the number of
     * slots available in the diagonal.  The sum along the diagonal will then be the count at that distance,
     * an "expected" or average uniform density.
     */
    public synchronized void computeDensity() {
	
		long maxNumBins = 0;
	
		//System.err.println("# of bins=" + numberOfBins);
		/**
		 * Genome wide binned possible distances
		 */
		double[] possibleDistances = new double[numberOfBins];
	
		for (Chromosome chr : chromosomesMap.values()) {
		
			// didn't see anything at all from a chromosome, then don't include it in possDists.
			if (chr == null || !chromosomeCounts.containsKey(chr.getIndex())) continue;
		
			// use correct units (bp or fragments)
            long len = chr.getLength();
            long nChrBins = len / gridSize;
		
			maxNumBins = Math.max(maxNumBins, nChrBins);
		
			for (int i = 0; i < nChrBins; i++) {
				possibleDistances[i] += (nChrBins - i);
			}
		
		}
		//System.err.println("max # bins " + maxNumBins);
		densityAvg = new ListOfDoubleArrays(maxNumBins);
	
		// Smoothing.  Keep pointers to window size.  When read counts drops below 400 (= 5% shot noise), smooth
        double shotNoiseMinimum = 400; //Math.log(1 + 400);
	
		double numSum = actualDistances[0];
		double denSum = possibleDistances[0];
		int bound1 = 0;
		int bound2 = 0;
		for (long ii = 0; ii < maxNumBins; ii++) {
            if (numSum < shotNoiseMinimum) {
                while (numSum < shotNoiseMinimum && bound2 < maxNumBins) {
                    // increase window size until window is big enough.  This code will only execute once;
                    // after this, the window will always contain at least 400 reads.
                    bound2++;
                    numSum += actualDistances[bound2];
                    denSum += possibleDistances[bound2];
                }
            } else if (numSum >= shotNoiseMinimum && bound2 - bound1 > 0) {
                while (bound2 - bound1 > 0 && bound2 < numberOfBins && bound1 < numberOfBins && numSum - actualDistances[bound1] - actualDistances[bound2] >= shotNoiseMinimum) {
                    numSum = numSum - actualDistances[bound1] - actualDistances[bound2];
                    denSum = denSum - possibleDistances[bound1] - possibleDistances[bound2];
                    bound1++;
                    bound2--;
                }
            }
            densityAvg.set(ii, numSum / denSum); // Math.expm1(numSum / denSum)
            // Default case - bump the window size up by 2 to keep it centered for the next iteration
            if (bound2 + 2 < maxNumBins) {
                numSum += actualDistances[bound2 + 1] + actualDistances[bound2 + 2];
                denSum += possibleDistances[bound2 + 1] + possibleDistances[bound2 + 2];
                bound2 += 2;
            } else if (bound2 + 1 < maxNumBins) {
                numSum += actualDistances[bound2 + 1];
                denSum += possibleDistances[bound2 + 1];
                bound2++;
            }
            // Otherwise, bound2 is at limit already
        }

        // Compute fudge factors for each chromosome so the total "expected" count for that chromosome == the observed

        for (Chromosome chr : chromosomesMap.values()) {
	
			if (chr == null || !chromosomeCounts.containsKey(chr.getIndex())) {
				continue;
			}
			//int len = isFrag ? fragmentCalculation.getNumberFragments(chr.getName()) : chr.getLength();
            long len = chr.getLength();
            long nChrBins = len / gridSize;
	
	
			double expectedCount = 0;
			for (long n = 0; n < nChrBins; n++) {
				if (n < maxNumBins) {
					final double v = densityAvg.get(n);
					// this is the sum of the diagonal for this particular chromosome.
					// the value in each bin is multiplied by the length of the diagonal to get expected count
					// the total at the end should be the sum of the expected matrix for this chromosome
					// i.e., for each chromosome, we calculate sum (genome-wide actual)/(genome-wide possible) == v
					// then multiply it by the chromosome-wide possible == nChrBins - n.
					expectedCount += (nChrBins - n) * v;
			
				}
            }

            double observedCount = chromosomeCounts.get(chr.getIndex());
            double f = expectedCount / observedCount;
            chrScaleFactors.put(chr.getIndex(), f);
        }
    }

    /**
     * Accessor for the normalization factors
     *
     * @return The normalization factors
     */
    public Map<Integer, Double> getChrScaleFactors() {
        return chrScaleFactors;
    }
	
	/**
	 * Accessor for the densities
	 *
	 * @return The densities
	 */
	public ListOfDoubleArrays getDensityAvg() {
		return densityAvg;
	}

    /**
     * Accessor for the normalization type
     *
     * @return The normalization type
     */
    public NormalizationType getType() {
        return type;
    }

    public ExpectedValueFunctionImpl getExpectedValueFunction() {
        computeDensity();
        return new ExpectedValueFunctionImpl(type, HiCZoom.HiCUnit.BP, gridSize, densityAvg, chrScaleFactors);
    }

    // TODO: this is often inefficient, we have all of the contact records when we leave norm calculations, should do this there if possible
    public void addDistancesFromIterator(int chrIndx, Iterator<ContactRecord> iterator, ListOfFloatArrays vector) {
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            int x = cr.getBinX();
            int y = cr.getBinY();
            final float counts = cr.getCounts();
            float xVal = vector.get(x);
            float yVal = vector.get(y);
            if (xVal > 0 & yVal > 0) {
                double value = counts / (xVal * yVal);
                addDistance(chrIndx, x, y, value);
            }
        }
    }
}


// Smooth in 3 stages,  the window sizes are tuned to human.

//        // Smooth (1)
//        final int smoothingWidow1 = 15000000;
//        int start = smoothingWidow1 / gridSize;
//        int window = (int) (5 * (2000000f / gridSize));
//        if (window == 0) window = 1;
//        for (int i = start; i < numberOfBins; i++) {
//            int kMin = i - window;
//            int kMax = Math.min(i + window, numberOfBins);
//            double sum = 0;
//            for (int k = kMin; k < kMax; k++) {
//                sum += density[k];
//            }
//            densityAvg[i] = sum / (kMax - kMin);
//        }
//
//        // Smooth (2)
//        start = 70000000 / gridSize;
//        window = (int)(20 * (2000000f / gridSize));
//        for (int i = start; i < numberOfBins; i++) {
//            int kMin = i - window;
//            int kMax = Math.min(i + window, numberOfBins);
//            double sum = 0;
//            for (int k = kMin; k < kMax; k++) {
//                sum += density[k];
//            }
//            densityAvg[i] = sum / (kMax - kMin);
//        }
//
//        // Smooth (3)
//        start = 170000000 / gridSize;
//        for (int i = start; i < numberOfBins; i++) {
//            densityAvg[i] = densityAvg[start];
//        }


/*

--- Code above based on the following Python

gridSize => grid (or bin) size  == 10^6
actualDistances => array of actual distances,  each element represents a bin
possibleDistances => array of possible distances, each element represents a bin
jdists => outer distances between pairs

for each jdist
  actualDistance[jdist]++;


for each chromosome
  chrlen = chromosome length
  numberOfBins = chrlen / gridSize
  for each i from 0 to numberOfBins
     possibleDistances[i] += (numberOfBins - i)


for each i from 0 to maxGrid
  density[i] = actualDistance[i] / possibleDistances[i]


for each i from 0 to len(density)
 density_avg[i] = density[i]

for each i from 15000000/gridsize  to  len(density_avg)
  sum1 = 0
  for each k from (i - 5*((2*10^6) / gridSize)  to  (i + 5*((2*10^6)/gridsize))
     sum1 += density[k]
  density_avg[i] = sum1 / (10*((2*10^6)/gridsize))

for each i from 70000000/gridsize  to  len(density_avg)
  sum2 = 0
  for each k from (i - 20*((2*10^6) / gridSize)  to  (i + 20*((2*10^6)/gridsize))
     sum2 += density[k]
  density_avg[i] = sum2 / (40*((2*10^6)/gridsize))

for each i from 170000000/gridsize  to  len(density_avg)
  density_avg[i]=density_avg[170000000/gridsize]

*/
