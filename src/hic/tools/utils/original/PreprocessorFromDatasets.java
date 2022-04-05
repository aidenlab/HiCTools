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

package hic.tools.utils.original;

import hic.tools.utils.iterators.contacts.AllByAllContactsIterator;
import hic.tools.utils.iterators.contacts.ChromosomeContactsIterator;
import hic.tools.utils.iterators.contacts.ContactIterator;
import hic.tools.utils.merge.HiCMergeTools;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.Deflater;


public class PreprocessorFromDatasets extends HiCFileBuilder {

    private static final Object key = new Object();

    public PreprocessorFromDatasets(File outputFile, String genomeId, ChromosomeHandler chromosomeHandler,
                                    double hicFileScalingFactor) {
        super(outputFile, genomeId, chromosomeHandler, hicFileScalingFactor);
    }

    public void preprocess(List<Dataset> inputDS, final String headerFile, final String footerFile) throws IOException {

        HiCMergeTools.mergeStatsAndGraphs(inputDS, tmpDir, this);

        try {
            LittleEndianOutputStream[] losFooter = initializeLosArrays(headerFile, footerFile);
            writeHeader();
            writeBody(inputDS);
            writeFooter(losFooter);
            closeLosArray(losFooter);
        } finally {
            closeLosArray(losArray);
        }

        updateMasterIndex(headerFile);
        System.out.println("\nFinished preprocess");
    }

    private MatrixPP computeWholeGenomeMatrix(List<Dataset> inputDS) {
        MatrixPP matrix = getInitialGenomeWideMatrixPP(chromosomeHandler);
        try {
            ContactIterator iter = new AllByAllContactsIterator(inputDS);
            while (iter.hasNext()) {
                matrix.incrementCount(iter.next(), expectedValueCalculations, tmpDir);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(98);
        }
        matrix.parsingComplete();
        return matrix;
    }

    private void writeBody(List<Dataset> inputDS) throws IOException {
        System.out.println("Writing body");
        writeWholeGenomeMatrix(inputDS, losArray, compressor, matrixPositions);

        Chromosome[] chromosomes = chromosomeHandler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                writeChromosomeRegionMatrix(chromosomes[i], chromosomes[j], inputDS);
            }
        }

        masterIndexPosition = losArray[0].getWrittenCount();
    }

    private void writeChromosomeRegionMatrix(Chromosome chromosome1, Chromosome chromosome2, List<Dataset> inputDS) throws IOException {
        int newBlockCapacity = (int) (Math.sqrt(inputDS.size()) * BLOCK_CAPACITY);
        final MatrixPP mergedMatrix = new MatrixPP(chromosome1.getIndex(), chromosome2.getIndex(), chromosomeHandler,
                bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity);

        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            MatrixPP matrixPP = new MatrixPP(chromosome1.getIndex(), chromosome2.getIndex(), chromosomeHandler,
                    bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity);

            while (i < inputDS.size()) {

                try {
                    ContactIterator iter = new ChromosomeContactsIterator(inputDS.get(i), chromosome1, chromosome2);
                    while (iter.hasNext()) {
                        matrixPP.incrementCount(iter.next(), expectedValueCalculations, tmpDir);
                    }
                } catch (Exception e) {
                    System.err.println("ERROR " + e.getLocalizedMessage());
                    System.err.println("Skipping dataset " + i + " for region: " +
                            chromosome1.getName() + "_" + chromosome2.getName());
                }

                i = index.getAndIncrement();
            }
            synchronized (key) {
                mergedMatrix.mergeMatrices(matrixPP);
            }
            matrixPP = null;
        });

        mergedMatrix.parsingComplete();
        writeMatrix(mergedMatrix, losArray, compressor, matrixPositions, -1, false);
    }

    private void writeWholeGenomeMatrix(List<Dataset> inputDS, LittleEndianOutputStream[] losArray, Deflater compressor,
                                        Map<String, IndexEntry> matrixPositions) throws IOException {
        MatrixPP wholeGenomeMatrix = computeWholeGenomeMatrix(inputDS);
        writeMatrix(wholeGenomeMatrix, losArray, compressor, matrixPositions, -1, false);
        wholeGenomeMatrix = null;
    }
}
