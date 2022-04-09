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
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.Deflater;


public class PreprocessorFromDatasets extends HiCFileBuilder {

    private static final Object key = new Object();
    private final Dataset[] datasets;
    private int highestResolution;

    public PreprocessorFromDatasets(File outputFile, Dataset[] datasets, double hicFileScalingFactor) {
        super(outputFile, datasets[0].getGenomeId(), hicFileScalingFactor);
        this.datasets = datasets;
        updateResolutionsToBuild(datasets[0].getAllPossibleResolutions());
        highestResolution = getMin(bpBinSizes);
    }

    private void updateResolutionsToBuild(List<HiCZoom> zooms) {
        List<Integer> resolutions = new ArrayList<>();
        for (HiCZoom zoom : zooms) {
            resolutions.add(zoom.getBinSize());
        }
        setResolutionsWithInts(resolutions);
    }

    private int getMin(int[] zooms) {
        int minValue = zooms[0];
        for (int zoom : zooms) {
            if (zoom < minValue) {
                minValue = zoom;
            }
        }
        return minValue;
    }

    public void preprocess() throws IOException {

        HiCMergeTools.mergeStatsAndGraphs(datasets, tmpDir, this);

        try {
            LittleEndianOutputStream[] losFooter = initializeLosArrays(outputFile.getAbsolutePath(),
                    outputFile.getAbsolutePath());
            writeHeader();
            writeBody();
            writeFooter(losFooter);
            closeLosArray(losFooter);
        } finally {
            closeLosArray(losArray);
        }

        updateMasterIndex(outputFile.getAbsolutePath());
        System.out.println("\nFinished preprocess");
    }

    private MatrixPP computeWholeGenomeMatrix(Dataset[] datasets) {
        MatrixPP matrix = getInitialGenomeWideMatrixPP(chromosomeHandler);
        try {
            ContactIterator iter = new AllByAllContactsIterator(datasets);
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

    private void writeBody() throws IOException {
        System.out.println("Writing body");
        writeWholeGenomeMatrix(datasets, losArray, compressor, matrixPositions);

        Chromosome[] chromosomes = chromosomeHandler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                if (diagonalsOnly && i != j) continue;
                writeChromosomeRegionMatrix(chromosomes[i], chromosomes[j], datasets);
                System.out.print(".");
            }
            System.out.println(".");
        }

        masterIndexPosition = losArray[0].getWrittenCount();
    }

    private void writeChromosomeRegionMatrix(Chromosome chromosome1, Chromosome chromosome2,
                                             Dataset[] datasets) throws IOException {
        int newBlockCapacity = (int) (Math.sqrt(datasets.length) * BLOCK_CAPACITY);
        final MatrixPP mergedMatrix = new MatrixPP(chromosome1.getIndex(), chromosome2.getIndex(), chromosomeHandler,
                bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity);

        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            MatrixPP matrixPP = new MatrixPP(chromosome1.getIndex(), chromosome2.getIndex(), chromosomeHandler,
                    bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity);

            while (i < datasets.length) {
                Matrix matrix = datasets[i].getMatrix(chromosome1, chromosome2, highestResolution);
                if (matrix == null) {
                    System.err.println("Skipping null matrix " + chromosome1.getName() + " " + chromosome2.getName());
                    continue;
                }
                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(highestResolution));
                if (zd == null) {
                    System.err.println("Skipping null zd (res=" + highestResolution + ") " +
                            chromosome1.getName() + " " + chromosome2.getName());
                    continue;
                }

                try {
                    ChromosomeContactsIterator iter = new ChromosomeContactsIterator(zd,
                            chromosome1, chromosome2, highestResolution);
                    if (!iter.hasNext()) {
                        System.err.println("No data in dataset " + i + " for region: " + zd.getKey());
                    }
                    while (iter.hasNext()) {
                        matrixPP.incrementCount(iter.next(), expectedValueCalculations, tmpDir);
                    }
                } catch (Exception e) {
                    System.err.println("ERROR " + e.getLocalizedMessage());
                    System.err.println("Skipping dataset " + i + " for region: " +
                            chromosome1.getName() + "_" + chromosome2.getName());
                    e.printStackTrace();
                    System.exit(90);
                }

                zd.clearCache();

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

    private void writeWholeGenomeMatrix(Dataset[] datasets, LittleEndianOutputStream[] losArray, Deflater compressor,
                                        Map<String, IndexEntry> matrixPositions) throws IOException {
        MatrixPP wholeGenomeMatrix = computeWholeGenomeMatrix(datasets);
        writeMatrix(wholeGenomeMatrix, losArray, compressor, matrixPositions, -1, false);
        wholeGenomeMatrix = null;
    }

    public void setHighestResolution(List<String> resolutions) {
        if (resolutions != null && resolutions.size() > 0) {
            try {
                highestResolution = Integer.parseInt(resolutions.get(0));
            } catch (Exception e) {
                System.err.println("Unable to parse resolution " + resolutions.get(0));
            }
        }
    }


}
