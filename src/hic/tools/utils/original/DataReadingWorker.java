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

package hic.tools.utils.original;

import hic.tools.utils.iterators.contacts.ChromosomeContactsIterator;
import hic.tools.utils.iterators.contacts.Contact;
import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;

import java.io.File;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class DataReadingWorker implements Runnable {

    private final AtomicInteger index;
    private final Chromosome chromosome1;
    private final Chromosome chromosome2;
    private final ChromosomeHandler chromosomeHandler;
    private final int[] bpBinSizes;
    private final int countThreshold;
    private final int v9DepthBase;
    private final int newBlockCapacity;
    private final int highestResolution;
    private final Dataset[] datasets;
    private final Map<String, ExpectedValueCalculation> expectedValueCalculations;
    private final File tmpDir;
    private final MatrixPP mergedMatrix;
    private final boolean onlyNearDiagonalContacts;

    public DataReadingWorker(AtomicInteger index, Chromosome chromosome1, Chromosome chromosome2,
                             ChromosomeHandler chromosomeHandler, int[] bpBinSizes, int countThreshold,
                             int v9DepthBase, int newBlockCapacity, Dataset[] datasets,
                             int highestResolution,
                             Map<String, ExpectedValueCalculation> expectedValueCalculations,
                             File tmpDir, MatrixPP mergedMatrix, boolean onlyNearDiagonalContacts) {
        this.index = index;
        this.chromosome1 = chromosome1;
        this.chromosome2 = chromosome2;
        this.chromosomeHandler = chromosomeHandler;
        this.bpBinSizes = bpBinSizes;
        this.countThreshold = countThreshold;
        this.v9DepthBase = v9DepthBase;
        this.newBlockCapacity = newBlockCapacity;
        this.datasets = datasets;
        this.highestResolution = highestResolution;
        this.expectedValueCalculations = expectedValueCalculations;
        this.tmpDir = tmpDir;
        this.mergedMatrix = mergedMatrix;
        this.onlyNearDiagonalContacts = onlyNearDiagonalContacts;
    }

    private static void processMatrix(Dataset[] datasets, int i, Chromosome chromosome1, Chromosome chromosome2,
                                      int highestResolution, boolean onlyNearDiagonalContacts, MatrixPP matrixPP,
                                      Map<String, ExpectedValueCalculation> expectedValueCalculations, File tmpDir) {
        Matrix matrix = datasets[i].getMatrix(chromosome1, chromosome2, highestResolution);
        if (matrix == null) {
            System.err.println("Skipping null matrix " + chromosome1.getName() + " " + chromosome2.getName());
            return;
        }
        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(highestResolution));
        if (zd == null) {
            System.err.println("Skipping null zd (res=" + highestResolution + ") " +
                    chromosome1.getName() + " " + chromosome2.getName());
            return;
        }

        try {
            ChromosomeContactsIterator iter = new ChromosomeContactsIterator(zd,
                    chromosome1, chromosome2, highestResolution);
            if (!iter.hasNext()) {
                System.err.println("No data in dataset " + i + " for region: " + zd.getKey());
            }
            while (iter.hasNext()) {
                Contact contact = iter.next();
                if (onlyNearDiagonalContacts && tooFarFromDiagonal(contact)) continue;
                matrixPP.incrementCount(contact, expectedValueCalculations, tmpDir);
            }
        } catch (Exception e) {
            System.err.println("ERROR " + e.getLocalizedMessage());
            System.err.println("Skipping dataset " + i + " for region: " +
                    chromosome1.getName() + "_" + chromosome2.getName());
            e.printStackTrace();
            System.exit(90);
        }

        zd.clearCache();
        datasets[i].clearCache(false);
    }

    private static boolean tooFarFromDiagonal(Contact contact) {
        return HiCFileBuilder.tooFarFromDiagonal(contact.getPos1(), contact.getPos2());
    }

    @Override
    public void run() {
        int i = index.getAndIncrement();
        MatrixPP matrixPP = new MatrixPP(chromosome1.getIndex(), chromosome2.getIndex(), chromosomeHandler,
                bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity);

        while (i < datasets.length) {
            processMatrix(datasets, i, chromosome1, chromosome2, highestResolution,
                    onlyNearDiagonalContacts, matrixPP, expectedValueCalculations, tmpDir);
            i = index.getAndIncrement();
        }
        synchronized (PreprocessorFromDatasets.key) {
            mergedMatrix.mergeMatrices(matrixPP);
        }
        matrixPP = null;
    }
}
