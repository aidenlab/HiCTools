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
import hic.tools.utils.largelists.BigShortsArray;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;

public class FinalScale {

    private final static short S1 = (short) 1;
    private final static short S0 = (short) 0;
    private final static float maxPercentile = 20f;
    private final static float startingPercentile = 1f;
    private final static float deltaPercentile = 0.5f;
    private final static float tolerance = 5e-4f;
    private final static int maxIter = 100;
    private final static int totalIterations = 3 * maxIter;
    private final static float minErrorThreshold = .02f;

    public static ListOfFloatArrays scaleToTargetVector(BigContactList ba, long matrixSize,
                                                        BigFloatsArray initialGuess, String stem) {

        long startTime = System.nanoTime();


        BigShortsArray bad = new BigShortsArray(matrixSize);

        BigShortsArray zTargetVector = new BigShortsArray(matrixSize, S1);
        BigFloatsArray calculatedVectorB = new BigFloatsArray(matrixSize);
        BigShortsArray one = new BigShortsArray(matrixSize, S1);

        double[] reportErrorForIteration = new double[totalIterations + 3];
        int[] numItersForAllIterations = new int[totalIterations + 3];

        //	find rows sums
        ListOfIntArrays numNonZero = ba.getNumNonZeroInRows();

        //	find relevant percentiles
        Percentiles percentiles = new Percentiles(numNonZero, startingPercentile, deltaPercentile, maxPercentile);
        int cutoffIndex = 0;
        int maxNumCutoffs = percentiles.getMaxNumCutoffs();
        double lowCutoff = percentiles.get(cutoffIndex);
        for (long p = 0; p < matrixSize; p++) {
            if (numNonZero.get(p) < lowCutoff) {
                excludeBadRow(p, bad, zTargetVector, one);
            }
        }

        BigFloatsArray row = ba.parSparseMultiplyAcrossLists(one, matrixSize);
        BigFloatsArray rowBackup = row.deepClone();

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
        //	start iterations
        //	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
        double convergeError = 10.0 * (1.0 + tolerance);
        double rowSumError = convergeError;
        int iter = 0;
        boolean failed;
        int nerr = 0;
        float[] errors = new float[10000];
        int allItersI = 0;
        BigFloatsArray col;
        int realIters = 0;

        // rowSumError = .1
        // total iter ~ 50

        // if perc or perc1 reached upper bound or the total number of iterations is too high, exit
        while ((convergeError > tolerance || rowSumError > 5.0 * tolerance)
                && iter < maxIter && allItersI < totalIterations
                && cutoffIndex < maxNumCutoffs) {
            iter++;
            allItersI++;
            realIters++;
            failed = true;

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

            reportErrorForIteration[allItersI - 1] = convergeError;
            numItersForAllIterations[allItersI - 1] = iter;

            //	since calculating the error in row sums requires matrix-vector multiplication
            //	we are doing this every 10	iterations
            if (iter % 10 == 0) {
                col = ba.parSparseMultiplyAcrossLists(calculatedVectorB, matrixSize);
                rowSumError = BigFloatsArray.parCalculateError(col, calculatedVectorB, zTargetVector, bad);
                errors[nerr++] = (float) rowSumError;
            }

            // current = calculatedVectorB
            current.parSetTo(calculatedVectorB);

            // check whether convergence rate is satisfactory
            // if less than 5 iterations (so less than 5 errors) and less than 2 row sums errors, there is nothing to check

            if (convergeError < tolerance
                    && (nerr < 2 || errors[nerr - 1] < 0.5 * errors[nerr - 2])) continue;

            if (iter > 5) {
                for (int q = 1; q <= 5; q++) {
                    if (reportErrorForIteration[allItersI - q] * (1.0 + minErrorThreshold) < reportErrorForIteration[allItersI - q - 1]) {
                        failed = false;
                        break;
                    }
                }

                if (nerr >= 2 && errors[nerr - 1] > 0.75 * errors[nerr - 2]) {
                    failed = true;
                }
    
                if (iter >= maxIter) {
                    failed = true;
                }

                if (failed) {
                    cutoffIndex++;
                    nerr = 0;
                    lowCutoff = percentiles.get(cutoffIndex);

                    for (long p = 0; p < matrixSize; p++) {
                        if (numNonZero.get(p) < lowCutoff && zTargetVector.get(p) > 0) {
                            excludeBadRow(p, bad, zTargetVector, one);
                        }
                    }

                    convergeError = 10.0 * (1.0 + tolerance);
                    rowSumError = 10.0 * (1.0 + tolerance);

                    //	if the current error is larger than 5 iteration ago start from scratch,
                    //	otherwise continue from the current position
                    if (reportErrorForIteration[allItersI - 1] > reportErrorForIteration[allItersI - 6]) {
                        for (long p = 0; p < matrixSize; p++) {
                            dr.set(p, 1 - bad.get(p));
                        }
                        dc.parSetTo(dr);
                        one.parSetTo(dr);
                        current.parSetTo(dr);
                        row.parSetTo(rowBackup);
                    } else {
                        // dr *= (1 - bad);
                        dr.parMultiplyByOneMinus(bad);
                        dc.parMultiplyByOneMinus(bad);
                    }
                    iter = 0;
                }
            }
        }

        //	find the final error in row sums
        col = ba.parSparseMultiplyAcrossLists(calculatedVectorB, matrixSize);
        rowSumError = BigFloatsArray.parCalculateError(col, calculatedVectorB, zTargetVector, bad);
        if (HiCGlobals.printVerboseComments) {
            double err90 = BigFloatsArray.calculateError90(col, calculatedVectorB, zTargetVector, bad);
            System.out.println("Total iters " + realIters + " \nRow Sums Error " + rowSumError + " \nError90 " + err90);
            System.out.println("Convergence error " + convergeError);
            reportErrorForIteration[allItersI + 1] = convergeError;
            reportErrorForIteration[allItersI + 2] = rowSumError;
        }

        if (convergeError > tolerance || rowSumError > 5.0 * tolerance || cutoffIndex >= maxNumCutoffs) {
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
        rowBackup.clear();
        dr.clear();
        dc.clear();
        current.clear();
        col.clear();

        if (HiCGlobals.printVerboseComments) {
            long endTime = System.nanoTime();
            long timeInSecs = (long) ((endTime - startTime) * 1e-9);

            System.out.println(stem + " took " + timeInSecs + " seconds");
            for (int q = 0; q < allItersI; q++) {
                System.out.println(numItersForAllIterations[q] + ": " + reportErrorForIteration[q]);
            }
            System.out.println("Total " + allItersI + " iterations; final zscore = " + percentiles.get(cutoffIndex));
            System.out.println("Final error in scaling vector is " + reportErrorForIteration[allItersI + 1] +
                    " and in row sums is " + reportErrorForIteration[allItersI + 2]);
        }

        ListOfFloatArrays answer = calculatedVectorB.convertToRegular();
        calculatedVectorB.clear();
        return answer;
    }

    private static void excludeBadRow(long index, BigShortsArray bad, BigShortsArray target, BigShortsArray one) {
        bad.set(index, S1);
        one.set(index, S0);
        target.set(index, S0);
    }

    private static BigFloatsArray update(long matrixSize, BigShortsArray bad,
                                         BigFloatsArray vector, BigShortsArray target,
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
