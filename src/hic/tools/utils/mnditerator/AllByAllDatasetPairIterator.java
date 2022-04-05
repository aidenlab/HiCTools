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

package hic.tools.utils.mnditerator;

import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * @author Jim Robinson
 * @since 4/7/12
 */
public class AllByAllDatasetPairIterator implements PairIterator {

    private final List<Dataset> datasets;
    private int currentIndex = -1;
    private Iterator<ContactRecord> iterator;
    private int resolution;

    public AllByAllDatasetPairIterator(List<Dataset> datasets) throws IOException {
        this.datasets = datasets;
        advanceIterator();
    }

    private static AlignmentPair makeAlignmentPairFromRecord(ContactRecord record, int resolution) {
        AlignmentPair pair = new AlignmentPair(0, resolution * record.getBinX(),
                0, resolution * record.getBinY());
        pair.setScore(record.getCounts());
        return pair;
    }

    private void advanceIterator() {
        currentIndex++;
        if (currentIndex < datasets.size()) {
            Chromosome all = datasets.get(currentIndex).getChromosomeHandler().getChromosomeFromIndex(0);
            Matrix matrix = datasets.get(currentIndex).getMatrix(all, all);
            MatrixZoomData zd = matrix.getFirstZoomData();
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
            if (currentIndex >= datasets.size()) {
                return false;
            }
        }
    }

    public AlignmentPair next() {
        ContactRecord record = iterator.next();
        return makeAlignmentPairFromRecord(record, resolution);
    }

    public void remove() {
    }

    public void close() {
    }
}
