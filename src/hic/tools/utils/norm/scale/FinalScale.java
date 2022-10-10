/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2020-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

package hic.tools.utils.norm.scale;

import hic.HiCGlobals;
import hic.tools.utils.bigarray.BigContactList;
import hic.tools.utils.largelists.BigFloatsArray;
import hic.tools.utils.largelists.BigIntsArray;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;

public class FinalScale {

    private final static short S1 = (short) 1;
    private final static short S0 = (short) 0;
    private final static float maxPercentile = 10f;
    private final static float tolerance = 1e-4f;
    private final static int maxIter = 500;
    private final static int totalIterations = 3 * maxIter;
    private final static float minErrorThreshold = .02f;
    private final static float del = .05f;
    private final static float rsError = .05f;
    private final static float erez = 1.0e-5f;


    public static ListOfFloatArrays scaleToTargetVector(BigContactList ba, long matrixSize,
                                                        BigFloatsArray initialGuess, String stem) {

        long startTime = System.nanoTime();

        int matrixSizeI = (int) matrixSize;
        BigIntsArray bad = new BigIntsArray(matrixSize);

        BigIntsArray zTargetVector = new BigIntsArray(matrixSize, S1);
        BigFloatsArray calculatedVectorB = new BigFloatsArray(matrixSize);
        BigIntsArray one = new BigIntsArray(matrixSize, S1);

        double[] reportErrorForIteration = new double[totalIterations + 3];
        int[] numItersForAllIterations = new int[totalIterations + 3];

        //	find numbers of nonzeros
        ListOfIntArrays numNonZero = ba.getNumNonZeroInRows();
        // find the highest number of nonzeros in a row we are willing to remove
        // may be better to use ListOfIntArrays instead, but it needs to be sorted
        int[] simpleNumNonZero = new int[matrixSizeI];
        int n0 = 0;
        for (long p = 0; p < matrixSize; p++) {
            int a = numNonZero.get(p);
            if (a > 0) simpleNumNonZero[n0++] = a;
        }
        java.util.Arrays.sort(simpleNumNonZero, 0, n0);
        // if we need to remove rows with more than upperBound nonzeros, scaling has failed
        int upperBound = simpleNumNonZero[(int) (maxPercentile * n0 / 100.0)];
//        delete(simpleNumNonZero);  // no longer needed

        int lowCutoff = 1;  // start with removing no nonzerow rows
        for (long p = 0; p < matrixSize; p++) {
            if (numNonZero.get(p) < lowCutoff) {
                excludeBadRow(p, bad, zTargetVector, one);
            }
        }

        BigFloatsArray row = ba.parSparseMultiplyAcrossLists(one, matrixSize);
    //    BigFloatsArray rowBackup = row.deepClone();

        for (long p = 0; p < matrixSize; p++) {
            one.set(p, (short) (1 - bad.get(p)));
        }

        BigFloatsArray dr;
        if (initialGuess == null) {
            dr = one.deepConvertedClone();
        } else {
            dr = initialGuess;
        }
        BigFloatsArray dc = dr.deepClone();
        BigFloatsArray current = dr.deepClone();
        row.parMultiplyBy(dr);

        // variables for the new algorithm
        boolean conv = false, div = false;
        int low_conv = 1000, low_div = 0;
        float[] b_conv = new float[matrixSizeI]; // keep the last converged vector
        float[] b0 = new float[matrixSizeI]; // keep the copy of the current vector
        int[] bad_conv = new int[matrixSizeI]; // bad rows for Erez's trick
        double ber_conv = 10.0;
        boolean yes = true;

        //	start iterations
        //	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
        double convergeError = 10.0 * (1.0 + tolerance);
        int iter = 0;
        int allItersI = 0;
        BigFloatsArray col;
        int realIters = 0;

        ITERS:
        while (convergeError > tolerance
                && iter < maxIter && allItersI < totalIterations) {
            iter++;
            allItersI++;
            realIters++;

            col = update(matrixSize, bad, row, zTargetVector, dr, ba);
            col.parMultiplyBy(dc);

            row = update(matrixSize, bad, col, zTargetVector, dc, ba);
            row.parMultiplyBy(dr);

            // calculate current scaling vector
            // calculatedVectorB[p] = (float) Math.sqrt(dr[p] * dc[p]);
            calculatedVectorB.parSetToGeoMean(dr, dc);
            //printFirst10(calculatedVectorB, iter);

            //	calculate the current error
            convergeError = BigFloatsArray.parCalculateConvergenceError(calculatedVectorB, current, bad);

            // new stuff
            double temp1;
            int numBad = 0;
            for (long p = 0; p < matrixSize; p++) {
                if (bad.get(p) == 1) continue;
                temp1 = Math.abs((calculatedVectorB.get(p) - current.get(p)) / (calculatedVectorB.get(p) + current.get(p)));
                if (temp1 > tolerance) numBad++;
            }
            for (int p = 0; p < matrixSizeI; p++) b0[p] = current.get(p);

            reportErrorForIteration[allItersI - 1] = convergeError;
            numItersForAllIterations[allItersI - 1] = iter;

            // current = calculatedVectorB
            current.parSetTo(calculatedVectorB);

            // check whether convergence rate is satisfactory
            // if less than 5 iterations (so less than 5 errors) and less than 2 row sums errors, there is nothing to check

            if (convergeError < tolerance) {  //if converged
                yes = true;
                if (lowCutoff == 1) break;
                conv = true;
                for (int p = 0; p < matrixSizeI; p++) b_conv[p] = calculatedVectorB.get(p);
                for (int p = 0; p < matrixSizeI; p++) bad_conv[p] = bad.get(p);
                ber_conv = convergeError;
                low_conv = lowCutoff;
                //  did it diverge before?
                if (div) {
                    if (low_conv - low_div <= 1) break;
                    lowCutoff = (low_conv + low_div) / 2;
                }
//  just halve low
                else lowCutoff = low_conv / 2;

                for (long p = 0; p < matrixSize; p++) {
                    one.set(p, S1);
                    bad.set(p, S0);
                }
                for (long p = 0; p < matrixSize; p++) {
                    if (numNonZero.get(p) < lowCutoff) {
                        bad.set(p, S1);
                        one.set(p, S0);
                    }
                }

                convergeError = 10.0;
                iter = 0;
                for (long p = 0; p < matrixSize; p++) dr.set(p, 1 - bad.get(p));
                for (long p = 0; p < matrixSize; p++) dc.set(p, 1 - bad.get(p));
                row = ba.parSparseMultiplyAcrossLists(dc, matrixSize);
                row.parMultiplyBy(dr);
                continue;
            }

            if (iter <= 5) continue;
//      check whether convergence rate is satisfactory
            if ((reportErrorForIteration[allItersI - 1] * (1.0 + del) < reportErrorForIteration[allItersI - 6]) && (iter < maxIter))
                continue;

//  diverged
            div = true;
            low_div = lowCutoff;
//  did it converge before? If it converged for low+1 and diverged for low, use the last converged norm vector
            if (conv) {
                if (low_conv - low_div <= 1) {
                    for (long p = 0; p < matrixSize; p++) calculatedVectorB.set(p, b_conv[(int) p]);
                    for (int p = 0; p < matrixSizeI; p++) bad.set(p, (short) bad_conv[p]);
                    convergeError = ber_conv;
                    break;
                }
//  if it almost converged (only a very small fraction of errors is above tol) remove bad rows and try again
//  with the same low (Erez's trick)
                else if (((double) numBad) / n0 < erez && yes) {
                    for (int p = 0; p < matrixSizeI; p++) {
                        if (bad.get(p) == 1) continue;
                        temp1 = Math.abs((calculatedVectorB.get(p) - b0[p]) / (calculatedVectorB.get(p) + b0[p]));
                        if (temp1 > tolerance) {
                            bad.set(p, S1);
                            one.set(p, S0);
                        }
                    }
                    yes = false;
                    convergeError = 10.0;
                    iter = 0;
                    for (long p = 0; p < matrixSize; p++) dr.set(p, 1.0f - bad.get(p));
                    for (long p = 0; p < matrixSize; p++) dc.set(p, 1.0f - bad.get(p));
                    row = ba.parSparseMultiplyAcrossLists(dc, matrixSize);
                    row.parMultiplyBy(dr);
                    //      if perc reached upper bound or the total number of iterationbs is too high, exit
                    if (lowCutoff > upperBound) break;
                    if (allItersI > totalIterations) break;
                    continue ITERS;
                } else {
                    lowCutoff = (low_div + low_conv) / 2;
                    yes = true;
                }
            }
            //  have never converged before
            //  Erez's trick
            else if (((double) numBad) / n0 < erez && yes) {
                for (int p = 0; p < matrixSizeI; p++) {
                    if (bad.get(p) == 1) continue;
                    temp1 = Math.abs((calculatedVectorB.get(p) - b0[p]) / (calculatedVectorB.get(p) + b0[p]));
                    if (temp1 > tolerance) {
                        bad.set(p, S1);
                        one.set(p, S0);
                    }
                }
                yes = false;
                convergeError = 10.0;
                iter = 0;
                for (long p = 0; p < matrixSize; p++) dr.set(p, 1.0f - bad.get(p));
                for (long p = 0; p < matrixSize; p++) dc.set(p, 1.0f - bad.get(p));
                row = ba.parSparseMultiplyAcrossLists(dc, matrixSize);
                row.parMultiplyBy(dr);
                //      if perc reached upper bound or the total number of iterationbs is too high, exit
                if (lowCutoff > upperBound) break;
                if (allItersI > totalIterations) break;
                continue ITERS;
            } else {
                lowCutoff = 2 * lowCutoff;
                yes = true;
            }

            for (long p = 0; p < matrixSize; p++) {
                bad.set(p, S0);
                one.set(p, S1);
            }
            for (long p = 0; p < matrixSize; p++) {
                if (numNonZero.get(p) < lowCutoff) {
                    bad.set(p, S1);
                    one.set(p, S0);
                }
            }
            convergeError = 10.0;
            iter = 0;
            for (long p = 0; p < matrixSize; p++) dr.set(p, 1.0f - bad.get(p));
            for (long p = 0; p < matrixSize; p++) dc.set(p, 1.0f - bad.get(p));
            row = ba.parSparseMultiplyAcrossLists(dc, matrixSize);
            row.parMultiplyBy(dr);
            //      if perc reached upper bound or the total number of iterationbs is too high, exit
            if (lowCutoff > upperBound) break;
            if (allItersI > totalIterations) break;
        }

        //	find the final error in row sums
        col = ba.parSparseMultiplyAcrossLists(calculatedVectorB, matrixSize);
        double rowSumError = BigFloatsArray.parCalculateError(col, calculatedVectorB, zTargetVector, bad);
        if (HiCGlobals.printVerboseComments) {
        //    double err90 = BigFloatsArray.calculateError90(col, calculatedVectorB, zTargetVector, bad);
        //    System.out.println("Total iters " + realIters + " \nRow Sums Error " + rowSumError + " \nError90 " + err90);
            System.out.println("Total iters " + realIters + " \nRow Sums Error " + rowSumError);
            System.out.println("Convergence error " + convergeError);
            System.out.println("Remove rows with less than " + lowCutoff + " nonzeros");
            reportErrorForIteration[allItersI + 1] = convergeError;
            reportErrorForIteration[allItersI + 2] = rowSumError;
        }

        if (convergeError > tolerance || rowSumError > rsError || lowCutoff > upperBound) {
            if (HiCGlobals.printVerboseComments) {
                System.out.println("Setting vector to null (not converged)");
            }
            calculatedVectorB.clear();
            return null;
        }

        for (long p = 0; p < matrixSize; p++) {
            if (bad.get(p) == 1) {
                calculatedVectorB.set(p, Float.NaN);
            }
        }

        bad.clear();
        zTargetVector.clear();
        one.clear();
        numNonZero.clear();
        row.clear();
        dr.clear();
        dc.clear();
        current.clear();
        col.clear();

        if (HiCGlobals.printVerboseComments) {
            long endTime = System.nanoTime();
            long timeInSecs = (long) ((endTime - startTime) * 1e-9);

            System.out.println(stem + " took " + timeInSecs + " seconds");
         /*
            for (int q = 0; q < allItersI; q++) {
                System.out.println(numItersForAllIterations[q] + ": " + reportErrorForIteration[q]);
            }
         */
            //System.out.println("Total " + allItersI + " iterations; final zscore = " + percentiles.get(cutoffIndex));
            System.out.println("Final error in scaling vector is " + reportErrorForIteration[allItersI + 1] +
                    " and in row sums is " + reportErrorForIteration[allItersI + 2]);
        }

        ListOfFloatArrays answer = calculatedVectorB.convertToRegular();
        calculatedVectorB.clear();
        return answer;
    }

    private static void excludeBadRow(long index, BigIntsArray bad, BigIntsArray target, BigIntsArray one) {
        bad.set(index, S1);
        one.set(index, S0);
        target.set(index, S0);
    }

    private static BigFloatsArray update(long matrixSize, BigIntsArray bad,
                                         BigFloatsArray vector, BigIntsArray target,
                                         BigFloatsArray dVector, BigContactList ba) {
        for (long p = 0; p < matrixSize; p++) {
            if (bad.get(p) == 1) {
                vector.set(p, 1f);
            }
        }
        // d *= target/vector
        dVector.parScaleByRatio(target, vector);

        // find sums and update scaling vector
        return ba.parSparseMultiplyAcrossLists(dVector, matrixSize);
    }
}
