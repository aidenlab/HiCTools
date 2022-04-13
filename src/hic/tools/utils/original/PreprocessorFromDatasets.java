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

import hic.tools.utils.iterators.contacts.AllByAllContactsIterator;
import hic.tools.utils.iterators.contacts.ChromosomeContactsIterator;
import hic.tools.utils.iterators.contacts.Contact;
import hic.tools.utils.iterators.contacts.ContactIterator;
import hic.tools.utils.merge.HiCMergeTools;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.Deflater;


public class PreprocessorFromDatasets extends HiCFileBuilder {

    private static final Object key = new Object();
    private final Dataset[] datasets;
    private int highestResolution;

    public PreprocessorFromDatasets(File outputFile, Dataset[] datasets, double hicFileScalingFactor) {
        super(outputFile, datasets[0].getGenomeId(), hicFileScalingFactor);
        this.datasets = datasets;
        updateResolutionsToBuild(datasets[0].getAllPossibleResolutions(), 0);
        highestResolution = getMin(bpBinSizes);
    }

    private void updateResolutionsToBuild(List<HiCZoom> zooms, int minCutoff) {
        List<Integer> resolutions = new ArrayList<>();
        for (HiCZoom zoom : zooms) {
            if (zoom.getBinSize() >= minCutoff) {
                resolutions.add(zoom.getBinSize());
            }
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

    private static boolean tooFarFromDiagonal(Contact contact) {
        return tooFarFromDiagonal(contact.getPos1(), contact.getPos2());
    }

    private Runnable getDataWritingWorker(LinkedList<MatrixPP> queue, LittleEndianOutputStream[] losArray, Deflater compressor,
                                          Map<String, IndexEntry> matrixPositions, AtomicBoolean stillHaveRegionsToProcess) {
        return () -> {
            while (queue.size() > 0 || stillHaveRegionsToProcess.get()) {
                if (queue.size() > 0) {
                    MatrixPP mergedMatrix = queue.pop();
                    try {
                        writeMatrix(mergedMatrix, losArray, compressor, matrixPositions,
                                -1, false);
                    } catch (IOException e) {
                        System.err.println("Unable to write matrix data to hic file");
                        e.printStackTrace();
                        System.exit(89);
                    }
                    mergedMatrix = null;
                }
            }
        };
    }

    private void writeBody() throws IOException {
        System.out.println("Writing body");
        writeWholeGenomeMatrix(datasets, losArray, compressor, matrixPositions);

        LinkedList<MatrixPP> queue = new LinkedList<>();
        AtomicBoolean stillHaveRegionsToProcess = new AtomicBoolean(true);

        ExecutorService writingExecutor = Executors.newFixedThreadPool(1);

        Runnable dataWritingWorker = getDataWritingWorker(queue, losArray, compressor, matrixPositions,
                stillHaveRegionsToProcess);
        writingExecutor.execute(dataWritingWorker);

        Chromosome[] chromosomes = chromosomeHandler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                if (intraChromosomalOnly && i != j) continue;
                tryToSleep(100);
                while (queue.size() > 2) {
                    tryToSleep(5000);
                }
                readInChromosomeRegionMatrix(chromosomes[i], chromosomes[j], datasets, queue);
                System.out.print("*" + queue.size() + "*");
            }
            System.out.println("*");
        }
        stillHaveRegionsToProcess.set(false);

        writingExecutor.shutdown();
        while (!writingExecutor.isTerminated()) {
        }

        masterIndexPosition = losArray[0].getWrittenCount();
    }

    private void tryToSleep(long time) {
        try {
            Thread.sleep(time);
        } catch (InterruptedException e) {
            System.err.println("Thread sleeping error " + e.getLocalizedMessage());
        }
    }

    private void readInChromosomeRegionMatrix(Chromosome chromosome1, Chromosome chromosome2,
                                              Dataset[] datasets, LinkedList<MatrixPP> queue) {
        int newBlockCapacity = (int) (Math.sqrt(datasets.length) * BLOCK_CAPACITY);
        final MatrixPP mergedMatrix = new MatrixPP(chromosome1.getIndex(), chromosome2.getIndex(), chromosomeHandler,
                bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity);

        AtomicInteger index = new AtomicInteger();

        int numCPUThreads = 4;
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);

        Runnable dataReadingWorker = getDataReadingWorker(index, chromosome1, chromosome2, chromosomeHandler,
                bpBinSizes, countThreshold, v9DepthBase, newBlockCapacity, datasets, highestResolution,
                expectedValueCalculations, tmpDir, mergedMatrix);

        for (int l = 0; l < numCPUThreads; ++l) {
            executor.execute(dataReadingWorker);
        }

        executor.shutdown();
        while (!executor.isTerminated()) {
        }

        mergedMatrix.parsingComplete();
        queue.add(mergedMatrix);
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
        updateResolutionsToBuild(datasets[0].getAllPossibleResolutions(), highestResolution);
    }

    private Runnable getDataReadingWorker(AtomicInteger index, Chromosome chromosome1, Chromosome chromosome2,
                                          ChromosomeHandler chromosomeHandler, int[] bpBinSizes, int countThreshold,
                                          int v9DepthBase, int newBlockCapacity, Dataset[] datasets,
                                          int highestResolution,
                                          Map<String, ExpectedValueCalculation> expectedValueCalculations,
                                          File tmpDir, MatrixPP mergedMatrix) {
        return () -> {
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

                i = index.getAndIncrement();
            }
            synchronized (key) {
                mergedMatrix.mergeMatrices(matrixPP);
            }
            matrixPP = null;
        };
    }
}
