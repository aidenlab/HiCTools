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

package juicebox.tools.utils.bigarray;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.iterators.GenomeWideIterator;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;

import java.util.Iterator;

public class BigArrayCreator {
    public static BigArray createFromZD(MatrixZoomData zd) {
        return populateBigArrayFromSingleIterator(zd.getDirectIterator(), 10000000, zd.getMatrixSize());
    }

    public static BigArray createForWholeGenome(Dataset dataset, ChromosomeHandler handler,
                                                HiCZoom zoom, boolean includeIntra) {
        int limit = 10000000;
        if (includeIntra) {
            limit = 50000000;
        }
        return populateBigArrayFromSingleIterator(new GenomeWideIterator(dataset, handler, zoom, includeIntra), limit,
                calculateGWSize(handler, zoom.getBinSize()));
    }

    public static BigArray populateBigArrayFromSingleIterator(Iterator<ContactRecord> iterator, int limit,
                                                              long matrixSize) {
        BigArray allRecords = new BigArray(limit, matrixSize);
        int[] x = new int[limit];
        int[] y = new int[limit];
        float[] c = new float[limit];
        int counter = 0;
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            x[counter] = cr.getBinX();
            y[counter] = cr.getBinY();
            c[counter] = cr.getCounts();
            counter++;
            if (counter > limit) {
                allRecords.addSubList(x, y, c);
                x = new int[limit];
                y = new int[limit];
                c = new float[limit];
                counter = 0;
            }
        }
        if (counter > 0) {
            allRecords.addSubList(x, y, c, counter);
        }
        return allRecords;
    }

    private static long calculateGWSize(ChromosomeHandler handler, int resolution) {
        long totalSize = 0;
        for (Chromosome c1 : handler.getChromosomeArrayWithoutAllByAll()) {
            totalSize += (c1.getLength() / resolution) + 1;
        }
        return totalSize;
    }

}
