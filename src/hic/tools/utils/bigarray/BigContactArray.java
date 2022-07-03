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

import com.google.common.util.concurrent.AtomicDouble;
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

    @Override
    public void clear() {
        binXs.clear();
        binYs.clear();
        binVals.clear();
    }

    private int getNumThreads() {
        return Math.min(HiCGlobals.normThreads, binXs.size());
    }

    @Override
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

    @Override
    public void clearIntraAndShiftInter() {
        return;
    }

    @Override
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

    @Override
    public ListOfFloatArrays getRowSums() {
        final ListOfFloatArrays totalRowSums = new ListOfFloatArrays(matrixSize, 0);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            ListOfFloatArrays sums = new ListOfFloatArrays(matrixSize);
            while (sIndx < binXs.size()) {
                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);
                float[] subBinVals = binVals.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    int x = subBinXs[z];
                    int y = subBinYs[z];
                    float value = subBinVals[z];
                    SparseMatrixTools.updateRowSums(sums, x, y, value);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (totalRowSums) {
                totalRowSums.addValuesFrom(sums);
            }
        });

        return totalRowSums;
    }

    @Override
    public double[] getNormMatrixSumFactor(ListOfFloatArrays norm) {
        final AtomicDouble matrixSum = new AtomicDouble(0);
        final AtomicDouble normSum = new AtomicDouble(0);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            double[] mSum = new double[1];
            double[] nSum = new double[1];
            while (sIndx < binXs.size()) {

                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);
                float[] subBinVals = binVals.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    int x = subBinXs[z];
                    int y = subBinYs[z];
                    float value = subBinVals[z];
                    SparseMatrixTools.sumScaleFactor(norm, mSum, nSum, x, y, value);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (matrixSum) {
                matrixSum.addAndGet(mSum[0]);
                normSum.addAndGet(nSum[0]);
            }
        });

        return new double[]{normSum.get(), matrixSum.get()};
    }

    public long getMatrixSize() {
        return matrixSize;
    }

    @Override
    public ListOfFloatArrays normalizeVectorByScaleFactor(ListOfFloatArrays newNormVector) {
        SparseMatrixTools.invertVector(newNormVector);

        final AtomicDouble normalizedSumTotal = new AtomicDouble(0);
        final AtomicDouble sumTotal = new AtomicDouble(0);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            double[] normSum = new double[1];
            double[] sum = new double[1];
            while (sIndx < binXs.size()) {

                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);
                float[] subBinVals = binVals.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    int x = subBinXs[z];
                    int y = subBinYs[z];
                    float counts = subBinVals[z];
                    SparseMatrixTools.sumRawAndNorm(normSum, sum, x, y, counts, newNormVector);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (normalizedSumTotal) {
                normalizedSumTotal.addAndGet(normSum[0]);
                sumTotal.addAndGet(sum[0]);
            }
        });

        double scaleFactor = Math.sqrt(normalizedSumTotal.get() / sumTotal.get());
        newNormVector.multiplyEverythingBy(scaleFactor);
        return newNormVector;
    }

    @Override
    public ListOfIntArrays getNumNonZeroInRows() {
        final ListOfIntArrays numNonZeros = new ListOfIntArrays(matrixSize);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            ListOfIntArrays nonZeros = new ListOfIntArrays(matrixSize);
            while (sIndx < binXs.size()) {

                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    int x = subBinXs[z];
                    int y = subBinYs[z];
                    nonZeros.addTo(x, 1);
                    if (x != y) {
                        nonZeros.addTo(y, 1);
                    }
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (numNonZeros) {
                numNonZeros.addValuesFrom(nonZeros);
            }
        });

        return numNonZeros;
    }

    @Override
    public void updateGenomeWideExpected(int chrIdx, ListOfFloatArrays vector, ExpectedValueCalculation exp) {
        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            float[] subBinVals = binVals.get(sIndx);

            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                float counts = subBinVals[z];
                SparseMatrixTools.populateNormedExpected(chrIdx, vector, exp, x, y, counts);
            }
        }
    }
}
