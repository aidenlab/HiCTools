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

import hic.HiCGlobals;
import hic.tools.utils.cleaner.ContactCleaner;
import hic.tools.utils.iterators.mnd.AlignmentPair;
import hic.tools.utils.iterators.mnd.AsciiPairIterator;
import hic.tools.utils.iterators.mnd.PairIterator;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.type.NormalizationHandler;
import javastraw.tools.ParallelizationTools;
import org.broad.igv.util.Pair;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.Deflater;


public class MultithreadedPreprocessor extends Preprocessor {
    public static final String CAT_SCRIPT = "_cat_outputs.sh";
    private final Map<Integer, String> chromosomePairIndexes = new ConcurrentHashMap<>();
    private final Map<String, Integer> chromosomePairIndexesReverse = new ConcurrentHashMap<>();
    private final Map<Integer, Integer> chromosomePairIndex1 = new ConcurrentHashMap<>();
    private final Map<Integer, Integer> chromosomePairIndex2 = new ConcurrentHashMap<>();
    private final int chromosomePairCounter;
    private final Map<Integer, Integer> nonemptyChromosomePairs = new ConcurrentHashMap<>();
    private final Map<Integer, Map<Integer, MatrixPP>> wholeGenomeMatrixParts = new ConcurrentHashMap<>();
    private final Map<String, IndexEntry> localMatrixPositions = new ConcurrentHashMap<>();
    private final Map<Integer, Long> matrixSizes = new ConcurrentHashMap<>();
    private final Map<Integer, Map<Long, List<IndexEntry>>> chromosomePairBlockIndexes;
    protected static int numCPUThreads = 1;
    private final Map<Integer, Map<String, ExpectedValueCalculation>> allLocalExpectedValueCalculations;
    protected static Map<Integer, List<Chunk>> mndIndex = null;
    private final AtomicInteger chunkCounter = new AtomicInteger(0);
    private int totalChunks = 0;
    private final ConcurrentHashMap<Integer, AtomicInteger> completedChunksPerChrPair = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, Integer> numChunksPerChrPair = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, AtomicInteger> chrPairCompleted = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, AtomicInteger> chrPairAvailableThreads = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, Integer> chrPairBlockCapacities = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, Integer> chunkCounterToChrPairMap = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, Integer> chunkCounterToChrChunkMap = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, Map<Integer, Pair<Pair<Integer, Integer>, MatrixPP>>> threadSpecificChrPairMatrices = new ConcurrentHashMap<>();
    private final ConcurrentHashMap<Integer, MatrixPP> finalChrMatrices = new ConcurrentHashMap<>();

    public MultithreadedPreprocessor(File outputFile, String genomeId, double hicFileScalingFactor, int numCPUThreads,
                                     String mndIndexFile, String tmpDir) throws IOException {
        super(outputFile, genomeId, hicFileScalingFactor, tmpDir);
        MultithreadedPreprocessor.numCPUThreads = numCPUThreads;
        chromosomeIndexes = MTIndexHandler.populateChromosomeIndexes(chromosomeHandler, numCPUThreads);
        chromosomePairCounter = MTIndexHandler.populateChromosomePairIndexes(chromosomeHandler,
                chromosomePairIndexes, chromosomePairIndexesReverse,
                chromosomePairIndex1, chromosomePairIndex2);
        setMndIndex(mndIndexFile, chromosomePairIndexes);
        this.chromosomePairBlockIndexes = new ConcurrentHashMap<>(chromosomePairCounter, (float) 0.75, numCPUThreads);
        this.allLocalExpectedValueCalculations = new ConcurrentHashMap<>(numCPUThreads, (float) 0.75, numCPUThreads);
    }

    public void setMndIndex(String mndIndexFile, Map<Integer, String> chromosomePairIndexes) throws IOException {
        if (mndIndexFile != null && mndIndexFile.length() > 1) {
            mndIndex = MTIndexHandler.readMndIndex(mndIndexFile, chromosomePairIndexes);
        } else {
            throw new IOException("No mndIndex provided");
        }
    }

    @Override
    public void preprocess(final String inputFile, String ignore1, String ignore2, Map<Integer,
            List<Chunk>> ignore3) throws IOException {
        super.preprocess(inputFile, outputFile + "_header", outputFile + "_footer", mndIndex);

        try {
            PrintWriter finalOutput = new PrintWriter(outputFile + CAT_SCRIPT);
            StringBuilder catOutputLine = new StringBuilder();
            StringBuilder removeLine = new StringBuilder();
            catOutputLine.append("cat ").append(outputFile).append("_header");
            removeLine.append("rm ").append(outputFile).append("_header");
            for (int i = 0; i < chromosomePairCounter; i++) {
                if ((nonemptyChromosomePairs.containsKey(i) && chromosomePairBlockIndexes.containsKey(i) && mndIndex.containsKey(i)) || i == 0) {
                    catOutputLine.append(" ").append(outputFile).append("_").append(chromosomePairIndexes.get(i));
                    removeLine.append(" ").append(outputFile).append("_").append(chromosomePairIndexes.get(i));
                    if (i > 0) {
                        int numOfNeededThreads = chrPairAvailableThreads.get(i).get();
                        if (numOfNeededThreads > 1) {
                            for (int j = 1; j <= numOfNeededThreads * numResolutions; j++) {
                                catOutputLine.append(" ").append(outputFile).append("_").append(chromosomePairIndexes.get(i)).append("_").append(j);
                                removeLine.append(" ").append(outputFile).append("_").append(chromosomePairIndexes.get(i)).append("_").append(j);
                            }
                        }
                    }
                }
            }
            catOutputLine.append(" ").append(outputFile).append("_footer").append(" > ").append(outputFile).append("\n");
            removeLine.append(" ").append(outputFile).append("_footer\n");
            finalOutput.println(catOutputLine);
            finalOutput.println(removeLine);
            finalOutput.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("Unable to write to catOutputs.sh");
            System.exit(70);
        }
    }

    private Pair<Pair<Integer,Integer>, MatrixPP> processIndividualMatrixChunk(String inputFile, int chunkNumber,
                                                                               int currentChrPair, Set<String> syncWrittenMatrices, Map<String, ExpectedValueCalculation>
                                                                                       localExpectedValueCalculations, int threadNum) throws IOException {
        MatrixPP wholeGenomeMatrix = getInitialGenomeWideMatrixPP(chromosomeHandler);
        int i = chunkNumber;
        int chunksProcessed = 0;

        String currentMatrixName;
        int currentPairIndex;

        int currentChr1 = -1;
        int currentChr2 = -1;
        MatrixPP currentMatrix = null;
        String currentMatrixKey;

        while (i < totalChunks) {
            int chrPair = chunkCounterToChrPairMap.get(i);
            if (chrPair != currentChrPair) {
                break;
            }
            int chrChunk = chunkCounterToChrChunkMap.get(i);
            List<Chunk> chunkPositions = mndIndex.get(chrPair);
            PairIterator iter = null;
            if (mndIndex == null) {
                System.err.println("No index for merged nodups file.");
                System.exit(67);
            } else {
                iter = new AsciiPairIterator(inputFile, chromosomeIndexes, chunkPositions.get(chrChunk),
                        chromosomeHandler);
            }
            ContactCleaner cleaner = new ContactCleaner(chromosomeHandler);

            while (iter.hasNext()) {
                AlignmentPair pair = iter.next();
                // skip pairs that mapped to contigs
                if (pair.isNotContigPair()) {
                    if (shouldSkipContact(pair)) continue;
                    // Flip pair if needed so chr1 < chr2
                    cleaner.updateLatestContact(pair);

                    // only increment if not intraFragment and passes the mapq threshold
                    if (cleaner.doesntMatchCurrentBlock(currentChr1, currentChr2)) {

                        // Start the next matrix
                        currentChr1 = cleaner.getChr1();
                        currentChr2 = cleaner.getChr2();
                        currentMatrixKey = currentChr1 + "_" + currentChr2;

                        currentMatrixName = cleaner.getCurrentMatrixName();

                        currentPairIndex = chromosomePairIndexesReverse.get(currentMatrixName);

                        if (currentPairIndex != currentChrPair) {
                            break;
                        }

                        if (syncWrittenMatrices.contains(currentMatrixKey)) {
                            System.err.println("Error: the chromosome combination " + currentMatrixKey + " appears in multiple blocks");
                            if (outputFile != null) outputFile.deleteOnExit();
                            System.exit(58);
                        }
                        currentMatrix = new MatrixPP(currentChr1, currentChr2, chromosomeHandler, bpBinSizes,
                                countThreshold, v9DepthBase, chrPairBlockCapacities.get(currentChrPair));
                    }
                    cleaner.incrementCount(currentMatrix, localExpectedValueCalculations, tmpDir);

                    cleaner.incrementGWCount(wholeGenomeMatrix, localExpectedValueCalculations, tmpDir);
                }
            }

            iter.close();
            chunksProcessed++;
            i = chunkCounter.getAndIncrement();
        }
        if (currentMatrix != null) {
            currentMatrix.parsingComplete();
        }
        wholeGenomeMatrixParts.get(currentChrPair).put(threadNum, wholeGenomeMatrix);
        return new Pair<>(new Pair<>(i, chunksProcessed), currentMatrix);

    }

    @Override
    protected void writeBody(String inputFile, Map<Integer, List<Chunk>> mndIndex) throws IOException {

        Set<String> syncWrittenMatrices = Collections.synchronizedSet(new HashSet<>());
        for (int chrPair = 1; chrPair < chromosomePairCounter; chrPair++) {
            if (mndIndex.containsKey(chrPair)) {
                int numOfChunks = mndIndex.get(chrPair).size();
                completedChunksPerChrPair.put(chrPair, new AtomicInteger(0));
                numChunksPerChrPair.put(chrPair, numOfChunks);
                chrPairCompleted.put(chrPair, new AtomicInteger(0));
                chrPairAvailableThreads.put(chrPair, new AtomicInteger(0));
                chrPairBlockCapacities.put(chrPair, BLOCK_CAPACITY/Math.min(numCPUThreads,numOfChunks));
                threadSpecificChrPairMatrices.put(chrPair, new ConcurrentHashMap<>());
                wholeGenomeMatrixParts.put(chrPair, new ConcurrentHashMap<>());
                for (int i=0; i<numOfChunks; i++) {
                    int currentChunk = totalChunks;
                    chunkCounterToChrPairMap.put(currentChunk, chrPair);
                    chunkCounterToChrChunkMap.put(currentChunk, i);
                    totalChunks++;
                }
            }
        }

        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int i = 1; i < numCPUThreads; i++) {
            int threadNum = i;
            Runnable worker = () -> {
                try {
                    int currentChunk = chunkCounter.getAndIncrement();
                    Map<String, ExpectedValueCalculation> localExpectedValueCalculations = new LinkedHashMap<>();
                    for (int bBinSize : bpBinSizes) {
                        ExpectedValueCalculation calc = new ExpectedValueCalculation(chromosomeHandler, bBinSize, NormalizationHandler.NONE);
                        String key = "BP_" + bBinSize;
                        localExpectedValueCalculations.put(key, calc);
                    }
                    while (currentChunk < totalChunks) {
                        int currentChrPair = chunkCounterToChrPairMap.get(currentChunk);
                        threadSpecificChrPairMatrices.get(currentChrPair).put(threadNum, processIndividualMatrixChunk(inputFile, currentChunk, currentChrPair, syncWrittenMatrices, localExpectedValueCalculations, threadNum));
                        synchronized (finalChrMatrices) {
                            if (!finalChrMatrices.containsKey(currentChrPair)) {
                                int currentChr1 = chromosomePairIndex1.get(currentChrPair);
                                int currentChr2 = chromosomePairIndex2.get(currentChrPair);
                                finalChrMatrices.put(currentChrPair, new MatrixPP(currentChr1, currentChr2, chromosomeHandler, bpBinSizes, countThreshold, v9DepthBase, chrPairBlockCapacities.get(currentChrPair)));
                            }
                            synchronized (finalChrMatrices.get(currentChrPair)) {
                                finalChrMatrices.get(currentChrPair).mergeMatrices(threadSpecificChrPairMatrices.get(currentChrPair).get(threadNum).getSecond());
                            }
                        }

                        for (int completedChunks = 0; completedChunks < threadSpecificChrPairMatrices.get(currentChrPair).get(threadNum).getFirst().getSecond(); completedChunks++) {
                            completedChunksPerChrPair.get(currentChrPair).getAndIncrement();
                        }
                        //System.err.println(currentChrPair + " " + threadSpecificChrPairMatrices.get(currentChrPair).get(threadNum).getFirst().getSecond() + " " + Duration.between(A,B).toMillis() + " " + Duration.between(B,C).toMillis() + " " + completedChunksPerChrPair.get(currentChrPair).get());
                        currentChunk = threadSpecificChrPairMatrices.get(currentChrPair).get(threadNum).getFirst().getFirst();
                        int currentAvailableThreads = chrPairAvailableThreads.get(currentChrPair).incrementAndGet();
                        if (completedChunksPerChrPair.get(currentChrPair).get() == numChunksPerChrPair.get(currentChrPair)) {
                            writeIndividualMatrix(currentChrPair, currentAvailableThreads);
                            finalChrMatrices.remove(currentChrPair);
                            threadSpecificChrPairMatrices.remove(currentChrPair);
                            chrPairCompleted.get(currentChrPair).getAndIncrement();
                            //System.err.println(currentChrPair + " " + Duration.between(D,E).toMillis());
                        }
                        while (chrPairCompleted.get(currentChrPair).get() == 0) {
                            try {
                                Thread.sleep(1000);
                            } catch (InterruptedException e) {
                                System.err.println(e.getLocalizedMessage());
                            }
                        }

                    }
                    allLocalExpectedValueCalculations.put(threadNum, localExpectedValueCalculations);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            };
            executor.execute(worker);
        }

        ParallelizationTools.shutDownAndWaitUntilDone(executor, 1000);

        for (int i = 0; i < numCPUThreads; i++) {
            if (allLocalExpectedValueCalculations.get(i) != null) {
                for (Map.Entry<String, ExpectedValueCalculation> entry : allLocalExpectedValueCalculations.get(i).entrySet()) {
                    expectedValueCalculations.get(entry.getKey()).merge(entry.getValue());
                }
            }
        }

        MatrixPP wholeGenomeMatrix = getInitialGenomeWideMatrixPP(chromosomeHandler);

        for (int i = 1; i < chromosomePairCounter; i++) {
            if (nonemptyChromosomePairs.containsKey(i)) {
                if (wholeGenomeMatrixParts.containsKey(i)) {
                    for (Map.Entry<Integer, MatrixPP> entry : wholeGenomeMatrixParts.get(i).entrySet()) {
                        wholeGenomeMatrix.mergeMatrices(entry.getValue());
                    }
                }
            }
        }

        // just making this more readable
        FileOutputStream tempFOS = new FileOutputStream(outputFile + "_" + chromosomePairIndexes.get(0));
        LittleEndianOutputStream tempLOS = new LittleEndianOutputStream(new BufferedOutputStream(tempFOS, HiCGlobals.bufferSize));
        LittleEndianOutputStream[] localLos = {tempLOS};
        writeMatrix2(wholeGenomeMatrix, localLos, WriterUtils.getDefaultCompressor(), localMatrixPositions,
                0, outputFile);
        nonemptyChromosomePairs.put(0, 1);

        long currentPosition = losArray[0].getWrittenCount();
        long nextMatrixPosition;
        String currentMatrixKey;

        for (int i = 0; i < chromosomePairCounter; i++) {
            if (nonemptyChromosomePairs.containsKey(i) && chromosomePairBlockIndexes.containsKey(i)) {
                for (Map.Entry<Long, List<IndexEntry>> entry : chromosomePairBlockIndexes.get(i).entrySet()) {
                    updateIndexPositions(entry.getValue(), null, false,
                            new File(outputFile + "_" + chromosomePairIndexes.get(i)),
                            currentPosition, entry.getKey());
                }
                nextMatrixPosition = localMatrixPositions.get("" + i).position + currentPosition;
                currentMatrixKey = chromosomePairIndex1.get(i) + "_" + chromosomePairIndex2.get(i);
                matrixPositions.put(currentMatrixKey, new IndexEntry(nextMatrixPosition, localMatrixPositions.get("" + i).size));
                currentPosition += matrixSizes.get(i);
            }
        }

        masterIndexPosition = currentPosition;
    }

    private void writeIndividualMatrix(Integer chromosomePair, int numOfNeededThreads) throws IOException {
        int chr1 = chromosomePairIndex1.get(chromosomePair);
        int chr2 = chromosomePairIndex2.get(chromosomePair);
        if (includedChromosomes != null) {
            String c1Name = chromosomeHandler.getChromosomeFromIndex(chr1).getName();
            String c2Name = chromosomeHandler.getChromosomeFromIndex(chr2).getName();
            if (includedChromosomes.contains(c1Name) || includedChromosomes.contains(c2Name)) {
                nonemptyChromosomePairs.put(chromosomePair, 1);
            }
        } else {
            nonemptyChromosomePairs.put(chromosomePair, 1);
        }

        LittleEndianOutputStream[] localLos;
        if (numOfNeededThreads == 1) {
            localLos = new LittleEndianOutputStream[1];
            localLos[0] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile + "_" + chromosomePairIndexes.get(chromosomePair)), HiCGlobals.bufferSize));
        } else {
            localLos = new LittleEndianOutputStream[(numOfNeededThreads * numResolutions) + 1];
            localLos[0] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile + "_" + chromosomePairIndexes.get(chromosomePair)), HiCGlobals.bufferSize));
            for (int i = 1; i <= numOfNeededThreads * numResolutions; i++) {
                localLos[i] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile + "_" + chromosomePairIndexes.get(chromosomePair) + "_" + i), HiCGlobals.bufferSize));
            }
        }

        writeMatrix2(finalChrMatrices.get(chromosomePair), localLos, WriterUtils.getDefaultCompressor(),
                localMatrixPositions, chromosomePair, outputFile);

    }

    // MatrixPP matrix, LittleEndianOutputStream los, Deflater compressor
    protected Pair<Map<Long, List<IndexEntry>>, Long> writeMatrix2(MatrixPP matrix, LittleEndianOutputStream[] localLos,
                                                                   Deflater localCompressor, Map<String, IndexEntry> localMatrixPositions,
                                                                   int chromosomePairIndex, File outputFile) throws IOException {

        Pair<Map<Long, List<IndexEntry>>, Long> localBlockIndexes = writeMatrix(matrix, localLos, localCompressor,
                localMatrixPositions, chromosomePairIndex, true, outputFile);

        chromosomePairBlockIndexes.put(chromosomePairIndex, localBlockIndexes.getFirst());
        long size = -localBlockIndexes.getSecond();
        for (LittleEndianOutputStream localLo : localLos) {
            size += localLo.getWrittenCount();
            localLo.close();
        }
        matrixSizes.put(chromosomePairIndex, size);
        return localBlockIndexes;
    }
}
