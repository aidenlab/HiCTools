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

package hic.tools.utils.bigarray;

import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.Iterator;

public class BigContactArrayCreator {
    public static BigContactArray createFromZD(MatrixZoomData zd) {
        return populateBigArrayFromSingleIterator(zd.getDirectIterator(), 10000000, zd.getMatrixSize());
    }

    public static BigContactArray populateBigArrayFromSingleIterator(Iterator<ContactRecord> iterator, int limit,
                                                                     long matrixSize) {
        BigContactArray allRecords = new BigContactArray(matrixSize);
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
            if (counter >= limit) {
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

}
