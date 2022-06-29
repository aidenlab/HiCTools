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

package hic.tools.utils.iterators.contacts;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;

import java.io.IOException;
import java.util.Iterator;

public class AllByAllContactsIterator implements ContactIterator {

    private final Dataset[] datasets;
    private int currentIndex = -1;
    private MatrixZoomData zd = null;

    private Matrix matrix = null;
    private Iterator<ContactRecord> iterator = null;
    private int resolution;

    public AllByAllContactsIterator(Dataset[] datasets) throws IOException {
        this.datasets = datasets;
        advanceIterator();
    }

    private static Contact makeAlignmentPairFromRecord(ContactRecord record, int resolution) {
        return new Contact(0, 0, record, resolution);
    }

    private void advanceIterator() {
        if (currentIndex > -1 && zd != null && matrix != null) {
            if (iterator != null) iterator = null;
            //matrix.clearCache();
            //datasets[currentIndex].clearCache(false);
        }

        currentIndex++;

        if (currentIndex < datasets.length) {
            Chromosome all = datasets[currentIndex].getChromosomeHandler().getChromosomeFromIndex(0);
            matrix = datasets[currentIndex].getMatrix(all, all);
            zd = matrix.getFirstZoomData();
            iterator = zd.getDirectIterator();
            resolution = zd.getBinSize();
        } else {
            iterator = null;
        }
    }

    public boolean hasNext() {
        if (iterator.hasNext()) {
            return true;
        }
        while (true) {
            advanceIterator();
            if (iterator != null && iterator.hasNext()) {
                return true;
            }
            if (currentIndex >= datasets.length) {
                return false;
            }
        }
    }

    public Contact next() {
        ContactRecord record = iterator.next();
        return makeAlignmentPairFromRecord(record, resolution);
    }
}
