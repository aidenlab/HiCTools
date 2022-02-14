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
import juicebox.tools.utils.bigarray.BigContactArray;
import juicebox.tools.utils.largelists.NormListOfFloatArrays;
import juicebox.tools.utils.largelists.NormListOfShortArrays;

import java.util.Arrays;

public class FinalScale {

    private final static short S1 = (short) 1;
    private final static short S0 = (short) 0;
    private final static float tol = .0005f;
    private final static float percentLowRowSumExcluded = 0.0001f;
    private final static float dp = percentLowRowSumExcluded / 2;
    private final static float tolerance = .0005f;
    private final static int maxIter = 100;
    private final static int totalIterations = 3 * maxIter;
    private final static float minErrorThreshold = .02f;
    private static final float OFFSET = .5f;

    public static ListOfFloatArrays scaleToTargetVector(BigContactArray ba, long matrixSize,
                                                        NormListOfFloatArrays initialGuess) {

        double low;
        int rlind;
        float localPercentLowRowSumExcluded = percentLowRowSumExcluded;
        NormListOfShortArrays bad = new NormListOfShortArrays(matrixSize);
        NormListOfFloatArrays s = new NormListOfFloatArrays(matrixSize);
        int[] r0 = new int[(int) Math.min(matrixSize, Integer.MAX_VALUE - 1)];

        NormListOfShortArrays zTargetVector = new NormListOfShortArrays(matrixSize, S1);
        NormListOfFloatArrays calculatedVectorB = new NormListOfFloatArrays(matrixSize);
        NormListOfShortArrays one = new NormListOfShortArrays(matrixSize, S1);

        double[] reportErrorForIteration = new double[totalIterations + 3];
        int[] numItersForAllIterations = new int[totalIterations + 3];

        //	find rows sums
        ListOfIntArrays numNonZero = ba.getNumNonZeroInRows();

        //	find relevant percentiles
        int n0 = 0;
        for (long p = 0; p < matrixSize; p++) {
            int valP = numNonZero.get(p);
            if (valP > 0 && n0 < r0.length) {
                r0[n0++] = valP;
            }
        }
        r0 = dealWithSorting(r0, n0);

        rlind = (int) Math.max(0, n0 * localPercentLowRowSumExcluded + OFFSET);
        low = r0[rlind];

        //	find the "bad" rows and exclude them
        for (long p = 0; p < matrixSize; p++) {
            if (numNonZero.get(p) < low) {
                bad.set(p, S1);
                zTargetVector.set(p, S0);
            }
        }

        NormListOfFloatArrays row = sparseMultiplyGetRowSums(ba, one, matrixSize);
        NormListOfFloatArrays rowBackup = row.deepClone();

        for (long p = 0; p < matrixSize; p++) {
            one.set(p, (short) (1 - bad.get(p)));
        }

        NormListOfFloatArrays dr;
        if (initialGuess == null) {
            dr = one.deepConvertedClone();
        } else {
            dr = initialGuess.deepClone();
        }
        NormListOfFloatArrays dc = dr.deepClone();
        NormListOfFloatArrays current = dr.deepClone();
        //	start iterations
        //	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
        double ber = 10.0 * (1.0 + tolerance);
        double err = ber;
        int iter = 0;
        boolean failed;
        int nerr = 0;
        double[] errors = new double[10000];
        int allItersI = 0;
        NormListOfFloatArrays col;

        // if perc or perc1 reached upper bound or the total number of iterationbs is too high, exit
        while ((ber > tolerance || err > 5.0 * tolerance) && iter < maxIter && allItersI < totalIterations
                && localPercentLowRowSumExcluded <= 0.2) {

            iter++;
            allItersI++;
            failed = true;

            col = update(matrixSize, bad, row, zTargetVector, s, dr, ba);
            multiply(col, dc, matrixSize);

            row = update(matrixSize, bad, col, zTargetVector, s, dc, ba);
            multiply(row, dr, matrixSize);

            // calculate current scaling vector
            for (long p = 0; p < matrixSize; p++) {
                calculatedVectorB.set(p, (float) Math.sqrt(dr.get(p) * dc.get(p)));
            }
            //printFirst10(calculatedVectorB, iter);

            //	calculate the current error
            ber = 0;
            for (long p = 0; p < matrixSize; p++) {
                if (bad.get(p) == 1) continue;
                double tempErr = Math.abs(calculatedVectorB.get(p) - current.get(p));
                if (tempErr > ber) {
                    ber = tempErr;
                }
            }
    
            reportErrorForIteration[allItersI - 1] = ber;
            numItersForAllIterations[allItersI - 1] = iter;

            //	since calculating the error in row sums requires matrix-vector multiplication we are doing this every 10
            //	iterations
            if (iter % 10 == 0) {
                col = sparseMultiplyGetRowSums(ba, calculatedVectorB, matrixSize);
                err = 0;
                for (long p = 0; p < matrixSize; p++) {
                    if (bad.get(p) == 1) continue;
                    double tempErr = Math.abs((col.get(p) * calculatedVectorB.get(p) - zTargetVector.get(p)));
                    if (err < tempErr) {
                        err = tempErr;
                    }
                }
                errors[nerr++] = err;
            }

            current.clear();
            current = calculatedVectorB.deepClone();

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
                    localPercentLowRowSumExcluded += dp;
                    nerr = 0;
                    rlind = (int) Math.max(0, n0 * localPercentLowRowSumExcluded + OFFSET);
                    low = r0[rlind];
                    for (long p = 0; p < matrixSize; p++) {
                        if (numNonZero.get(p) < low && zTargetVector.get(p) > 0) {
                            bad.set(p, S1);
                            one.set(p, S0);
                            zTargetVector.set(p, S0);
                        }
                    }

                    ber = 10.0 * (1.0 + tol);
                    err = 10.0 * (1.0 + tol);
    
                    //	if the current error is larger than 5 iteration ago start from scratch,
                    //	otherwise continue from the current position
                    if (reportErrorForIteration[allItersI - 1] > reportErrorForIteration[allItersI - 6]) {
                        if (initialGuess != null) {
                            for (long p = 0; p < matrixSize; p++) {
                                dr.set(p, (1 - bad.get(p)) * initialGuess.get(p));
                            }
                        } else {
                            // todo dr.setOneMinus(bad);
                            for (long p = 0; p < matrixSize; p++) {
                                dr.set(p, 1 - bad.get(p));
                            }
                        }
                        dc.clear();
                        one.clear();
                        current.clear();
                        row.clear();
                        dc = dr.deepClone();
                        one = dr.deepCovertedClone();
                        current = dr.deepClone();
                        row = rowBackup.deepClone();
                    } else {
                        for (long p = 0; p < matrixSize; p++) {
                            dr.multiplyBy(p, (1 - bad.get(p)));
                        }
                        for (long p = 0; p < matrixSize; p++) {
                            dc.multiplyBy(p, (1 - bad.get(p)));
                        }
                    }
                    iter = 0;
                }
            }
        }

        //	find the final error in row sums
        /*
        if (iter % 10 == 0) {
            col = sparseMultiplyGetRowSums(ba, calculatedVectorB, matrixSize);
            err = 0;
            for (int p = 0; p < matrixSize; p++) {
                if (bad.get(p) == 1) continue;
                double tempErr = Math.abs(col.get(p) * calculatedVectorB.get(p) - zTargetVector.get(p));
                if (err < tempErr) {
                    err = tempErr;
                }
            }
        }
        
        
        reportErrorForIteration[allItersI + 1] = ber;
        reportErrorForIteration[allItersI + 2] = err;
        */

        for (long p = 0; p < matrixSize; p++) {
            if (bad.get(p) == 1) {
                calculatedVectorB.set(p, Float.NaN);
            }
        }

        /*
        if (HiCGlobals.printVerboseComments) {
            System.out.println(allItersI);
            System.out.println(localPercentLowRowSumExcluded);
            System.out.println(Arrays.toString(reportErrorForIteration));
        }
        */

        //printFirst10(calculatedVectorB, -1);
        //System.out.println("DONE");
        return calculatedVectorB.convertToRegular();
    }

    private static void printFirst10(ListOfFloatArrays calculatedVectorB, int i) {
        float[] row = new float[10];
        for (int z = 0; z < 10; z++) {
            row[z] = calculatedVectorB.get(z);
        }
        System.out.println("ITER " + i + " NORM VEC " + Arrays.toString(row));
    }

    private static NormListOfFloatArrays update(long matrixSize, NormListOfShortArrays bad,
                                                NormListOfFloatArrays vector, NormListOfShortArrays target,
                                                NormListOfFloatArrays s, NormListOfFloatArrays dVector,
                                                BigContactArray ic) {
        for (long p = 0; p < matrixSize; p++) {
            if (bad.get(p) == 1) vector.set(p, 1.0f);
        }
        for (long p = 0; p < matrixSize; p++) {
            s.set(p, target.get(p) / vector.get(p));
        }
        multiply(dVector, s, matrixSize);

        // find sums and update scaling vector
        return sparseMultiplyGetRowSums(ic, dVector, matrixSize);
    }

    private static void multiply(NormListOfFloatArrays vec, NormListOfFloatArrays dv, long matrixSize) {
        for (long p = 0; p < matrixSize; p++) {
            vec.multiplyBy(p, dv.get(p));
        }
    }

    private static int[] dealWithSorting(int[] vector, int length) {
        int[] realVector = new int[length];
        System.arraycopy(vector, 0, realVector, 0, length);
        Arrays.sort(realVector);
        return realVector;
    }

    private static NormListOfFloatArrays sparseMultiplyGetRowSums(BigContactArray ba, NormListOfFloatArrays vector, long vectorLength) {
        return ba.sparseMultiplyAcrossLists(vector, vectorLength);
    }

    private static NormListOfFloatArrays sparseMultiplyGetRowSums(BigContactArray ba, NormListOfShortArrays vector, long vectorLength) {
        return ba.sparseMultiplyAcrossLists(vector, vectorLength);
    }
}
