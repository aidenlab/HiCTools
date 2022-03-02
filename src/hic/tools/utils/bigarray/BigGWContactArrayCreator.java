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

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.iterators.GenomeWideIterator;
import javastraw.reader.type.HiCZoom;

import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class BigGWContactArrayCreator {
    public static BigGWContactArray createForWholeGenome(Dataset dataset, ChromosomeHandler handler,
                                                         HiCZoom zoom, boolean includeIntra) {
        int limit = 10000000;
        if (includeIntra) {
            limit = 50000000;
        }
        long matrixSize = calculateGWSize(handler, zoom.getBinSize());

        BigContactArray[] bas = parLoadInterAndIntraGWData(dataset, handler, zoom, includeIntra, limit, matrixSize);

        BigGWContactArray array = new BigGWContactArray(matrixSize);
        array.addAllSubLists(bas[0]);
        array.addAllSubListsToInterBackup(bas[0]);
        bas[0].clear();
        if (bas[1] != null) {
            array.addAllSubLists(bas[1]);
            bas[1].clear();
        }
        bas = null;

        return array;
    }

    private static BigContactArray[] parLoadInterAndIntraGWData(Dataset dataset, ChromosomeHandler handler, HiCZoom zoom, boolean includeIntra, int limit, long matrixSize) {
        BigContactArray[] bas = new BigContactArray[2];
        Arrays.fill(bas, null);
        ExecutorService executor = Executors.newFixedThreadPool(2);

        Runnable worker = () -> {
            bas[0] = BigContactArrayCreator.populateBigArrayFromSingleIterator(new GenomeWideIterator(dataset, handler,
                    zoom, false, true), limit, matrixSize);
        };
        executor.execute(worker);

        Runnable worker2 = () -> {
            if (includeIntra) {
                bas[1] = BigContactArrayCreator.populateBigArrayFromSingleIterator(new GenomeWideIterator(dataset, handler,
                        zoom, true, false), limit, matrixSize);
            }
        };
        executor.execute(worker2);

        executor.shutdown();

        while (!executor.isTerminated()) {
        }

        return bas;
    }

    private static long calculateGWSize(ChromosomeHandler handler, int resolution) {
        long totalSize = 0;
        for (Chromosome c1 : handler.getChromosomeArrayWithoutAllByAll()) {
            totalSize += (c1.getLength() / resolution) + 1;
        }
        return totalSize;
    }
}
