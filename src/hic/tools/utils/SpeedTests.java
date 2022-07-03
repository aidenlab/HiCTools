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

package hic.tools.utils;

import hic.HiCGlobals;
import hic.tools.utils.bigarray.BigContactArrayCreator;
import hic.tools.utils.bigarray.BigContactList;
import javastraw.reader.Dataset;
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;

import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;

public class SpeedTests {

    public static void testRowSums() throws IOException {

        String path = "/Users/muhammad/Desktop/hicfiles/copy_tmp_chr10_subsample0.25.hic";

        DatasetReaderV2 reader = new DatasetReaderV2(path, false, false);
        Dataset ds = reader.read();

        Chromosome chr10 = ds.getChromosomeHandler().getChromosomeFromName("chr10");

        Matrix matrix = ds.getMatrix(chr10, chr10, 50);
        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(50));

        HiCGlobals.normThreads = 10;
        double[] timeTotals = new double[3];

        for (int i = 0; i < 4; i++) {
            // test 1 direct iteration
            ListOfFloatArrays f0 = test0(zd, timeTotals);
            System.gc();

            // test 1 big contact iterator
            ListOfFloatArrays f1 = test1(zd, timeTotals);
            System.gc();

            // test 1 local file based iterator
            ListOfFloatArrays f2 = test2(zd, timeTotals);
            System.gc();

            assertAreEqual(f0, f1, "f0_vs_f1");
            assertAreEqual(f0, f2, "f0_vs_f2");
        }

        System.out.println(Arrays.toString(timeTotals));

    }

    private static ListOfFloatArrays test0(MatrixZoomData zd, double[] timeTotals) {

        long r0 = System.nanoTime();
        ListOfFloatArrays f1 = testOnIterator(zd);
        long r1 = System.nanoTime();
        double time = (r1 - r0) * 1e-9;
        timeTotals[0] += time;
        System.err.println("\nTest 0: " + time + " seconds");
        return f1;
    }

    private static ListOfFloatArrays testOnIterator(MatrixZoomData zd) {
        final ListOfFloatArrays sums = new ListOfFloatArrays(zd.getMatrixSize());

        Iterator<ContactRecord> iterator = zd.getDirectIterator();
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            int x = record.getBinX();
            int y = record.getBinY();
            float value = record.getCounts();
            sums.addTo(x, value);
            if (x != y) {
                sums.addTo(y, value);
            }
        }

        return sums;
    }

    private static ListOfFloatArrays test1(MatrixZoomData zd, double[] timeTotals) {
        BigContactList ba = BigContactArrayCreator.createFromZD(zd);
        long r0 = System.nanoTime();
        ListOfFloatArrays f1 = ba.getRowSums();
        long r1 = System.nanoTime();
        double time = (r1 - r0) * 1e-9;
        timeTotals[1] += time;
        System.err.println("\nTest 1: " + time + " seconds");
        return f1;
    }

    private static ListOfFloatArrays test2(MatrixZoomData zd, double[] timeTotals) {
        BigContactList ba = BigContactArrayCreator.createLocalVersionFromZD(zd);
        long r0 = System.nanoTime();
        ListOfFloatArrays f1 = ba.getRowSums();
        long r1 = System.nanoTime();
        double time = (r1 - r0) * 1e-9;
        timeTotals[2] += time;
        System.err.println("\nTest 2: " + time + " seconds");
        return f1;
    }

    public static void assertAreEqual(ListOfFloatArrays data1, ListOfFloatArrays data2, String description) {
        double magnitude = 0;
        double absError = 0;
        try {
            if (data1.getLength() != data2.getLength()) {
                System.err.println("Vector length mismatch: " + data1.getLength() + " vs " + data2.getLength() + " " + description);
                //System.exit(24);
            }

            long n = Math.min(data1.getLength(), data2.getLength());

            double magnitude1 = 0;
            double magnitude2 = 0;

            long numVals = 0;
            for (long q = 0; q < n; q++) {
                double err = Math.abs(data1.get(q) - data2.get(q));

                if (!Double.isNaN(data1.get(q))) {
                    magnitude1 += data1.get(q) * data1.get(q);
                }

                if (!Double.isNaN(data2.get(q))) {
                    magnitude2 += data2.get(q) * data2.get(q);
                }

                if (!Double.isNaN(err)) {
                    magnitude += data1.get(q) * data2.get(q);
                    absError += err;
                    numVals++;
                }
            }
            magnitude = Math.sqrt(magnitude);
            magnitude1 = Math.sqrt(magnitude1);
            magnitude2 = Math.sqrt(magnitude2);

            if (numVals > 0) {
                absError /= numVals;
            }

            System.err.println("Vector mean error: (" + absError + ") / " + magnitude + " / " + magnitude1
                    + " / " + magnitude2 + " / " + description);

        } catch (Exception e) {
            System.exit(26);
        }
    }
}
