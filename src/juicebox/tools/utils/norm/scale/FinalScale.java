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

package juicebox.tools.utils.norm.scale;

import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;
import juicebox.HiCGlobals;
import juicebox.tools.utils.bigarray.BigContactArray;
import juicebox.tools.utils.largelists.BigFloatsArray;
import juicebox.tools.utils.largelists.BigShortsArray;

import java.util.Arrays;

public class FinalScale {

    private final static short S1 = (short) 1;
    private final static short S0 = (short) 0;
    private final static float tol = .0005f;
    private final static float zscoreCutoff = -3;
    private final static float dZscore = 0.25f;
    private final static float tolerance = .0005f;
    private final static int maxIter = 100;
    private final static int totalIterations = 3 * maxIter;
    private final static float minErrorThreshold = .02f;

    public static ListOfFloatArrays scaleToTargetVector(BigContactArray ba, long matrixSize,
                                                        BigFloatsArray initialGuess) {

        float localZscoreCutoff = zscoreCutoff;
        BigShortsArray bad = new BigShortsArray(matrixSize);
        BigFloatsArray svec = new BigFloatsArray(matrixSize);
        int[] r0 = new int[(int) Math.min(matrixSize, Integer.MAX_VALUE - 1)];

        BigShortsArray zTargetVector = new BigShortsArray(matrixSize, S1);
        BigFloatsArray calculatedVectorB = new BigFloatsArray(matrixSize);
        BigShortsArray one = new BigShortsArray(matrixSize, S1);

        double[] reportErrorForIteration = new double[totalIterations + 3];
        int[] numItersForAllIterations = new int[totalIterations + 3];

        //	find rows sums
        ListOfIntArrays numNonZero = ba.getNumNonZeroInRows();

        //	find relevant percentiles
        ZScore zscore = new ZScore(numNonZero);
        float lowCutoff = zscore.getCutoff(localZscoreCutoff);
        for (long p = 0; p < matrixSize; p++) {
            if (numNonZero.get(p) < lowCutoff) {
                excludeBadRow(p, bad, zTargetVector, one);
            }
        }

        BigFloatsArray row = parSparseMultiplyGetRowSums(ba, one, matrixSize);
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
        double ber = 10.0 * (1.0 + tolerance);
        double err = ber;
        int iter = 0;
        boolean failed;
        int nerr = 0;
        float[] errors = new float[10000];
        int allItersI = 0;
        BigFloatsArray col = null;
        int realIters = 0;

        // if perc or perc1 reached upper bound or the total number of iterationbs is too high, exit
        while ((ber > tolerance || err > 5.0 * tolerance) && iter < maxIter && allItersI < totalIterations
                && localZscoreCutoff <= -.5) {
            iter++;
            allItersI++;
            realIters++;
            failed = true;

            col = update(matrixSize, bad, row, zTargetVector, svec, dr, ba);
            col.parMultiplyBy(dc);

            row = update(matrixSize, bad, col, zTargetVector, svec, dc, ba);
            row.parMultiplyBy(dr);

            // calculate current scaling vector
            // calculatedVectorB[p] = (float) Math.sqrt(dr[p] * dc[p]);
            calculatedVectorB.parSetToGeoMean(dr, dc);
            //printFirst10(calculatedVectorB, iter);

            //	calculate the current error
            ber = BigFloatsArray.parCalculateConvergenceError(calculatedVectorB, current, bad);


            reportErrorForIteration[allItersI - 1] = ber;
            numItersForAllIterations[allItersI - 1] = iter;

            //	since calculating the error in row sums requires matrix-vector multiplication we are doing this every 10
            //	iterations
            if (iter % 10 == 0) {
                col = parSparseMultiplyGetRowSums(ba, calculatedVectorB, matrixSize);
                err = BigFloatsArray.parCalculateError(col, calculatedVectorB, zTargetVector, bad);
                errors[nerr++] = (float) err;
            }

            // current = calculatedVectorB
            current.parSetTo(calculatedVectorB);

            // check whether convergence rate is satisfactory
            // if less than 5 iterations (so less than 5 errors) and less than 2 row sums errors, there is nothing to check

            if (ber < tolerance && (nerr < 2 || errors[nerr - 1] < 0.5 * errors[nerr - 2])) continue;

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
                    localZscoreCutoff += dZscore;
                    nerr = 0;
                    lowCutoff = zscore.getCutoff(localZscoreCutoff);
                    for (long p = 0; p < matrixSize; p++) {
                        if (numNonZero.get(p) < lowCutoff && zTargetVector.get(p) > 0) {
                            excludeBadRow(p, bad, zTargetVector, one);
                        }
                    }

                    ber = 10.0 * (1.0 + tol);
                    err = 10.0 * (1.0 + tol);
    
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
        if (HiCGlobals.printVerboseComments) {
            col = parSparseMultiplyGetRowSums(ba, calculatedVectorB, matrixSize);
            err = BigFloatsArray.parCalculateError(col, calculatedVectorB, zTargetVector, bad);
            double err90 = BigFloatsArray.calculateError90(col, calculatedVectorB, zTargetVector, bad);
            System.out.println("Total iters " + realIters + " Row Sums Error " + err + " Error90 " + err90);

            reportErrorForIteration[allItersI + 1] = ber;
            reportErrorForIteration[allItersI + 2] = err;

            //System.out.println(allItersI);
            //System.out.println(localPercentLowRowSumExcluded);
            //System.out.println(Arrays.toString(reportErrorForIteration));
            //printFirst10(calculatedVectorB, -1);
        }

        for (long p = 0; p < matrixSize; p++) {
            if (bad.get(p) == 1) {
                calculatedVectorB.set(p, Float.NaN);
            }
        }

        bad.clear();
        svec.clear();
        zTargetVector.clear();
        one.clear();
        numNonZero.clear();
        r0 = null;
        row.clear();
        rowBackup.clear();
        dr.clear();
        dc.clear();
        current.clear();
        if (col != null) {
            col.clear();
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

    private static void printFirst10(ListOfFloatArrays calculatedVectorB, int i) {
        float[] row = new float[10];
        for (int z = 0; z < 10; z++) {
            row[z] = calculatedVectorB.get(z);
        }
        System.out.println("ITER " + i + " NORM VEC " + Arrays.toString(row));
    }

    private static BigFloatsArray update(long matrixSize, BigShortsArray bad,
                                         BigFloatsArray vector, BigShortsArray target,
                                         BigFloatsArray s, BigFloatsArray dVector,
                                         BigContactArray ic) {
        for (long p = 0; p < matrixSize; p++) {
            if (bad.get(p) == 1) {
                vector.set(p, 1f);
            }
        }
        // s = target/vector
        s.parSetToDivision(target, vector);
        dVector.parMultiplyBy(s);

        // find sums and update scaling vector
        return parSparseMultiplyGetRowSums(ic, dVector, matrixSize);
    }

    private static int[] dealWithSorting(int[] vector, int length) {
        int[] realVector = new int[length];
        System.arraycopy(vector, 0, realVector, 0, length);
        Arrays.sort(realVector);
        return realVector;
    }

    private static BigFloatsArray parSparseMultiplyGetRowSums(BigContactArray ba, BigFloatsArray vector, long vectorLength) {
        return ba.parSparseMultiplyAcrossLists(vector, vectorLength);
    }

    private static BigFloatsArray parSparseMultiplyGetRowSums(BigContactArray ba, BigShortsArray vector, long vectorLength) {
        return ba.parSparseMultiplyAcrossLists(vector, vectorLength);
    }
}
