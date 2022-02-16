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
import javastraw.reader.iterators.GenomeWideIterator;
import javastraw.reader.type.HiCZoom;

public class BigGWContactArrayCreator {
    public static BigGWContactArray createForWholeGenome(Dataset dataset, ChromosomeHandler handler,
                                                         HiCZoom zoom, boolean includeIntra) {
        int limit = 10000000;
        if (includeIntra) {
            limit = 50000000;
        }
        long matrixSize = calculateGWSize(handler, zoom.getBinSize());

        BigGWContactArray array = new BigGWContactArray(matrixSize);

        BigContactArray interBA = BigContactArrayCreator.populateBigArrayFromSingleIterator(new GenomeWideIterator(dataset, handler,
                        zoom, false, true), limit,
                matrixSize);
        array.addAllSubLists(interBA);
        array.addAllSubListsToInterBackup(interBA);
        interBA.clear();

        if (includeIntra) {
            BigContactArray intraBA = BigContactArrayCreator.populateBigArrayFromSingleIterator(new GenomeWideIterator(dataset, handler,
                    zoom, true, false), limit, matrixSize);
            array.addAllSubLists(intraBA);
            intraBA.clear();
        }

        return array;
    }

    private static long calculateGWSize(ChromosomeHandler handler, int resolution) {
        long totalSize = 0;
        for (Chromosome c1 : handler.getChromosomeArrayWithoutAllByAll()) {
            totalSize += (c1.getLength() / resolution) + 1;
        }
        return totalSize;
    }
}
