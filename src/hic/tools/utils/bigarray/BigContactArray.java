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

package hic.tools.utils.bigarray;

import hic.HiCGlobals;
import hic.tools.utils.largelists.BigDoublesArray;
import hic.tools.utils.largelists.BigFloatsArray;
import hic.tools.utils.largelists.BigShortsArray;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;
import javastraw.tools.ParallelizationTools;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class BigContactArray implements BigContactList {

    protected final List<int[]> binXs = new ArrayList<>();
    protected final List<int[]> binYs = new ArrayList<>();
    protected final List<float[]> binVals = new ArrayList<>();
    private final long matrixSize;

    public BigContactArray(long matrixSize) {
        this.matrixSize = matrixSize;
    }

    public void addSubList(int[] x, int[] y, float[] c) {
        binXs.add(x);
        binYs.add(y);
        binVals.add(c);
    }

    public void addSubList(int[] x, int[] y, float[] c, int counter) {
        int[] x2 = new int[counter];
        int[] y2 = new int[counter];
        float[] c2 = new float[counter];
        System.arraycopy(x, 0, x2, 0, counter);
        System.arraycopy(y, 0, y2, 0, counter);
        System.arraycopy(c, 0, c2, 0, counter);
        addSubList(x2, y2, c2);
    }

    public void addAllSubLists(BigContactArray other) {
        binXs.addAll(other.binXs);
        binYs.addAll(other.binYs);
        binVals.addAll(other.binVals);
    }

    public void clear() {
        binXs.clear();
        binYs.clear();
        binVals.clear();
    }

    private int getNumThreads() {
        return Math.min(HiCGlobals.normThreads, binXs.size());
    }

    public BigFloatsArray parSparseMultiplyAcrossLists(BigFloatsArray vector, long vectorLength) {
        final BigDoublesArray totalSumVector = new BigDoublesArray(vectorLength);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            BigDoublesArray sumVector = new BigDoublesArray(vectorLength);
            while (sIndx < binXs.size()) {
                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);
                float[] subBinVals = binVals.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    SparseMatrixTools.matrixVectorMult(vector, sumVector,
                            subBinXs[z], subBinYs[z], subBinVals[z]);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (totalSumVector) {
                totalSumVector.addValuesFrom(sumVector);
            }
        });

        return totalSumVector.convertToFloats();
    }

    public BigFloatsArray parSparseMultiplyAcrossLists(BigShortsArray vector, long vectorLength) {
        final BigDoublesArray totalSumVector = new BigDoublesArray(vectorLength);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            BigDoublesArray sumVector = new BigDoublesArray(vectorLength);
            while (sIndx < binXs.size()) {
                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);
                float[] subBinVals = binVals.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    SparseMatrixTools.matrixVectorMult(vector, sumVector,
                            subBinXs[z], subBinYs[z], subBinVals[z]);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (totalSumVector) {
                totalSumVector.addValuesFrom(sumVector);
            }
        });

        return totalSumVector.convertToFloats();
    }

    public ListOfFloatArrays getRowSums() {
        ListOfFloatArrays rowsums = new ListOfFloatArrays(matrixSize, 0);

        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            float[] subBinVals = binVals.get(sIndx);

            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                float value = subBinVals[z];
                rowsums.addTo(x, value);
                if (x != y) {
                    rowsums.addTo(y, value);
                }
            }
        }

        return rowsums;
    }

    public double[] getNormMatrixSumFactor(ListOfFloatArrays norm) {
        double matrix_sum = 0;
        double norm_sum = 0;
        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            float[] subBinVals = binVals.get(sIndx);

            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                float value = subBinVals[z];
                double valX = norm.get(x);
                double valY = norm.get(y);
                if (valX > 0 && valY > 0) {
                    // want total sum of matrix, not just upper triangle
                    if (x == y) {
                        norm_sum += value / (valX * valY);
                        matrix_sum += value;
                    } else {
                        norm_sum += 2 * value / (valX * valY);
                        matrix_sum += 2 * value;
                    }
                }
            }
        }
        return new double[]{norm_sum, matrix_sum};
    }

    public long getMatrixSize() {
        return matrixSize;
    }

    public ListOfFloatArrays normalizeVectorByScaleFactor(ListOfFloatArrays newNormVector) {
        SparseMatrixTools.invertVector(newNormVector);

        double normalizedSumTotal = 0, sumTotal = 0;

        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            float[] subBinVals = binVals.get(sIndx);

            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                float counts = subBinVals[z];

                double valX = newNormVector.get(x);
                double valY = newNormVector.get(y);

                if (valX > 0 && valY > 0) {
                    double normalizedValue = counts / (valX * valY);
                    normalizedSumTotal += normalizedValue;
                    sumTotal += counts;
                    if (x != y) {
                        normalizedSumTotal += normalizedValue;
                        sumTotal += counts;
                    }
                }
            }
        }

        double scaleFactor = Math.sqrt(normalizedSumTotal / sumTotal);
        newNormVector.multiplyEverythingBy(scaleFactor);
        return newNormVector;
    }

    public ListOfIntArrays getNumNonZeroInRows() {
        ListOfIntArrays numNonZero = new ListOfIntArrays(matrixSize, 0);
        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                numNonZero.addTo(x, 1);
                if (x != y) {
                    numNonZero.addTo(y, 1);
                }
            }
        }
        return numNonZero;
    }

    public void updateGenomeWideExpected(int chrIdx, ListOfFloatArrays vector, ExpectedValueCalculation exp) {
        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            float[] subBinVals = binVals.get(sIndx);

            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                float counts = subBinVals[z];
                if (counts > 0) {
                    final float vx = vector.get(x);
                    final float vy = vector.get(y);
                    if (vx > 0 && vy > 0) {
                        double value = counts / (vx * vy);
                        if (value > 0) {
                            exp.addDistance(chrIdx, x, y, value);
                        }
                    }
                }
            }
        }
    }
}
