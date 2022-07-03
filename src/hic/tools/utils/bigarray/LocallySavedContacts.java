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
import hic.tools.utils.localtemps.BinRecordsReader;
import hic.tools.utils.localtemps.BinRecordsWriter;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class LocallySavedContacts implements BigContactList {
    private final List<String> filenames = new ArrayList<>();
    private final long matrixSize;

    public LocallySavedContacts(Iterator<ContactRecord> directIterator, int bufferSize, long matrixSize) {
        this.matrixSize = matrixSize;
        try {
            BinRecordsWriter.saveAllContacts(directIterator, bufferSize, filenames);
        } catch (Exception e) {
            System.err.println("ERROR: Unable to save data locally");
            e.printStackTrace();
            System.exit(20);
        }
    }


    public void clear() {
        for (String filename : filenames) {
            File file = new File(filename);
            file.delete();
        }
        filenames.clear();
    }

    public long getMatrixSize() {
        return matrixSize;
    }

    private int getNumThreads() {
        return Math.min(HiCGlobals.normThreads, filenames.size());
    }

    public BigFloatsArray parSparseMultiplyAcrossLists(BigFloatsArray vector, long vectorLength) {
        final BigDoublesArray totalSumVector = new BigDoublesArray(vectorLength);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            BigDoublesArray sumVector = new BigDoublesArray(vectorLength);
            while (sIndx < filenames.size()) {

                try {
                    BinRecordsReader reader = new BinRecordsReader(filenames.get(sIndx));
                    while (reader.hasNext()) {
                        ContactRecord record = reader.next();
                        SparseMatrixTools.matrixVectorMult(vector, sumVector,
                                record.getBinX(), record.getBinY(), record.getCounts());
                    }
                    reader.close();
                } catch (IOException e) {
                    System.err.println("ERROR: unable to read file (" + filenames.get(sIndx) + "): " + e.getMessage());
                    e.printStackTrace();
                    System.exit(9);
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
            while (sIndx < filenames.size()) {

                try {
                    BinRecordsReader reader = new BinRecordsReader(filenames.get(sIndx));
                    while (reader.hasNext()) {
                        ContactRecord record = reader.next();
                        SparseMatrixTools.matrixVectorMult(vector, sumVector,
                                record.getBinX(), record.getBinY(), record.getCounts());

                    }
                    reader.close();
                } catch (IOException e) {
                    System.err.println("ERROR: unable to read file (" + filenames.get(sIndx) + "): " + e.getMessage());
                    e.printStackTrace();
                    System.exit(10);
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
        final ListOfFloatArrays totalRowSums = new ListOfFloatArrays(matrixSize, 0);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            ListOfFloatArrays sums = new ListOfFloatArrays(matrixSize);
            while (sIndx < filenames.size()) {
                try {
                    BinRecordsReader reader = new BinRecordsReader(filenames.get(sIndx));
                    while (reader.hasNext()) {
                        ContactRecord record = reader.next();
                        SparseMatrixTools.updateRowSums(sums, record.getBinX(), record.getBinY(), record.getCounts());
                    }
                    reader.close();
                } catch (IOException e) {
                    System.err.println("ERROR: unable to read file (" + filenames.get(sIndx) + "): " + e.getMessage());
                    e.printStackTrace();
                    System.exit(11);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (totalRowSums) {
                totalRowSums.addValuesFrom(sums);
            }
        });

        return totalRowSums;
    }

    public double[] getNormMatrixSumFactor(ListOfFloatArrays norm) {
        final AtomicDouble matrixSum = new AtomicDouble(0);
        final AtomicDouble normSum = new AtomicDouble(0);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            double[] mSum = new double[1];
            double[] nSum = new double[1];
            while (sIndx < filenames.size()) {

                try {
                    BinRecordsReader reader = new BinRecordsReader(filenames.get(sIndx));
                    while (reader.hasNext()) {
                        ContactRecord record = reader.next();

                        int x = record.getBinX();
                        int y = record.getBinY();
                        float value = record.getCounts();
                        SparseMatrixTools.sumScaleFactor(norm, mSum, nSum, x, y, value);
                    }
                    reader.close();
                } catch (IOException e) {
                    System.err.println("ERROR: unable to read file (" + filenames.get(sIndx) + "): " + e.getMessage());
                    e.printStackTrace();
                    System.exit(12);
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


    public ListOfFloatArrays normalizeVectorByScaleFactor(ListOfFloatArrays newNormVector) {
        SparseMatrixTools.invertVector(newNormVector);

        final AtomicDouble normalizedSumTotal = new AtomicDouble(0);
        final AtomicDouble sumTotal = new AtomicDouble(0);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            double[] normSum = new double[1];
            double[] sum = new double[1];
            while (sIndx < filenames.size()) {

                try {
                    BinRecordsReader reader = new BinRecordsReader(filenames.get(sIndx));
                    while (reader.hasNext()) {
                        ContactRecord record = reader.next();
                        int x = record.getBinX();
                        int y = record.getBinY();
                        float counts = record.getCounts();
                        SparseMatrixTools.sumRawAndNorm(normSum, sum, x, y, counts, newNormVector);
                    }
                    reader.close();
                } catch (IOException e) {
                    System.err.println("ERROR: unable to read file (" + filenames.get(sIndx) + "): " + e.getMessage());
                    e.printStackTrace();
                    System.exit(13);
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


    public ListOfIntArrays getNumNonZeroInRows() {
        final ListOfIntArrays numNonZeros = new ListOfIntArrays(matrixSize);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int sIndx = index.getAndIncrement();
            ListOfIntArrays nonZeros = new ListOfIntArrays(matrixSize);
            while (sIndx < filenames.size()) {
                try {
                    BinRecordsReader reader = new BinRecordsReader(filenames.get(sIndx));
                    while (reader.hasNext()) {
                        ContactRecord record = reader.next();
                        nonZeros.addTo(record.getBinX(), 1);
                        if (record.getBinX() != record.getBinY()) {
                            nonZeros.addTo(record.getBinY(), 1);
                        }
                    }
                    reader.close();
                } catch (IOException e) {
                    System.err.println("ERROR: unable to read file (" + filenames.get(sIndx) + "): " + e.getMessage());
                    e.printStackTrace();
                    System.exit(14);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (numNonZeros) {
                numNonZeros.addValuesFrom(nonZeros);
            }
        });

        return numNonZeros;
    }

    public void updateGenomeWideExpected(int chrIdx, ListOfFloatArrays expectedVector, ExpectedValueCalculation exp) {

        for (String filename : filenames) {
            try {
                BinRecordsReader reader = new BinRecordsReader(filename);
                while (reader.hasNext()) {
                    ContactRecord record = reader.next();
                    int x = record.getBinX();
                    int y = record.getBinY();
                    float counts = record.getCounts();
                    SparseMatrixTools.populateNormedExpected(chrIdx, expectedVector, exp, x, y, counts);
                }
                reader.close();
            } catch (IOException e) {
                System.err.println("ERROR: unable to read file (" + filename + "): " + e.getMessage());
                e.printStackTrace();
                System.exit(15);
            }
        }
    }
}
