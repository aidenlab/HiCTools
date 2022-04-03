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

//import juicebox.MainWindow;

import hic.HiCGlobals;
import hic.tools.utils.cleaner.ContactCleaner;
import hic.tools.utils.largelists.BigListOfByteWriters;
import hic.tools.utils.mnditerator.AlignmentPair;
import hic.tools.utils.mnditerator.PairIterator;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.Pair;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.Deflater;


/**
 * @author jrobinso
 * @since Aug 16, 2010
 */
public class Preprocessor extends HiCFileBuilder {

    /**
     * The position of the field containing the masterIndex position
     */
    protected long masterIndexPositionPosition;
    protected long normVectorIndexPosition;
    protected long normVectorLengthPosition;
    protected Map<String, ExpectedValueCalculation> expectedValueCalculations;

    public Preprocessor(File outputFile, String genomeId, ChromosomeHandler chromosomeHandler, double hicFileScalingFactor) {
        super(outputFile, genomeId, chromosomeHandler, hicFileScalingFactor);
    }

    public void preprocess(final String inputFile, final String headerFile, final String footerFile,
                           Map<Integer, List<Chunk>> mndIndex) throws IOException {
        File file = new File(inputFile);
        if (!file.exists() || file.length() == 0) {
            System.err.println(inputFile + " does not exist or does not contain any reads.");
            System.exit(57);
        }

        try {
            StringBuilder stats = null;
            StringBuilder graphs = null;
            StringBuilder hicFileScaling = new StringBuilder().append(hicFileScalingFactor);

            if (statsFileName != null) {
                try (FileInputStream is = new FileInputStream(statsFileName)) {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is), HiCGlobals.bufferSize);
                    stats = new StringBuilder();
                    String nextLine;
                    while ((nextLine = reader.readLine()) != null) {
                        stats.append(nextLine).append("\n");
                    }
                } catch (IOException e) {
                    System.err.println("Error while reading stats file: " + e);
                    stats = null;
                }

            }
            if (graphFileName != null) {
                try (FileInputStream is = new FileInputStream(graphFileName)) {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is), HiCGlobals.bufferSize);
                    graphs = new StringBuilder();
                    String nextLine;
                    while ((nextLine = reader.readLine()) != null) {
                        graphs.append(nextLine).append("\n");
                    }
                } catch (IOException e) {
                    System.err.println("Error while reading graphs file: " + e);
                    graphs = null;
                }
            }

            expectedValueCalculations = Collections.synchronizedMap(new LinkedHashMap<>());
            for (int bBinSize : bpBinSizes) {
                ExpectedValueCalculation calc = new ExpectedValueCalculation(chromosomeHandler, bBinSize, NormalizationHandler.NONE);
                String key = "BP_" + bBinSize;
                expectedValueCalculations.put(key, calc);
            }

            LittleEndianOutputStream[] losFooter = new LittleEndianOutputStream[1];
            try {
                losArray[0] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(headerFile), HiCGlobals.bufferSize));
                if (footerFile.equalsIgnoreCase(headerFile)) {
                    losFooter = losArray;
                } else {
                    losFooter[0] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(footerFile), HiCGlobals.bufferSize));
                }
            } catch (Exception e) {
                System.err.println("Unable to write to " + outputFile);
                System.exit(70);
            }

            System.out.println("Start preprocess");

            System.out.println("Writing header");

            writeHeader(stats, graphs, hicFileScaling);

            System.out.println("Writing body");
            writeBody(inputFile, mndIndex);

            System.out.println();
            System.out.println("Writing footer");
            writeFooter(losFooter);

            if (losFooter != null && losFooter[0] != null) {
                losFooter[0].close();
            }

        } finally {
            if (losArray != null && losArray[0] != null) {
                losArray[0].close();
            }
        }

        updateMasterIndex(headerFile);
        System.out.println("\nFinished preprocess");
    }

    protected void writeHeader(StringBuilder stats, StringBuilder graphs, StringBuilder hicFileScaling) throws IOException {
        // Magic number
        byte[] magicBytes = "HIC".getBytes();
        LittleEndianOutputStream los = losArray[0];
        los.write(magicBytes[0]);
        los.write(magicBytes[1]);
        los.write(magicBytes[2]);
        los.write(0);

        // VERSION
        los.writeInt(VERSION);
    
        // Placeholder for master index position, replaced with actual position after all contents are written
        masterIndexPositionPosition = los.getWrittenCount();
        los.writeLong(0L);
    
    
        // Genome ID
        los.writeString(genomeId);
    
        // Add NVI info
        //los.writeString(NVI_INDEX);
        normVectorIndexPosition = los.getWrittenCount();
        los.writeLong(0L);
    
        //los.writeString(NVI_LENGTH);
        normVectorLengthPosition = los.getWrittenCount();
        los.writeLong(0L);
    
    
        // Attribute dictionary
        int nAttributes = 1;
        if (stats != null) nAttributes += 1;
        if (graphs != null) nAttributes += 1;
        if (hicFileScaling != null) nAttributes += 1;
        if (v9DepthBase != 2) nAttributes += 1;

        los.writeInt(nAttributes);
        los.writeString(Dataset.SOFTWARE);
        los.writeString("Juicer Tools Version " + HiCGlobals.versionNum);
        if (stats != null) {
            los.writeString(Dataset.STATISTICS);
            los.writeString(stats.toString());
        }
        if (graphs != null) {
            los.writeString(Dataset.GRAPHS);
            los.writeString(graphs.toString());
        }
        if (hicFileScaling != null) {
            los.writeString(Dataset.HIC_FILE_SCALING);
            los.writeString(hicFileScaling.toString());
        }
        if (v9DepthBase != 2) {
            los.writeString(Dataset.V9_DEPTH_BASE);
            los.writeString("" + v9DepthBase);
        }


        // Sequence dictionary
        int nChrs = chromosomeHandler.size();
        los.writeInt(nChrs);
        for (Chromosome chromosome : chromosomeHandler.getChromosomeArray()) {
            los.writeString(chromosome.getName());
            los.writeLong(chromosome.getLength());
        }

        //BP resolution levels
        int nBpRes = bpBinSizes.length;
        los.writeInt(nBpRes);
        for (int bpBinSize : bpBinSizes) {
            los.writeInt(bpBinSize);
        }

        // deprecate fragment resolutions
        los.writeInt(0);

        numResolutions = nBpRes;
    }

    protected MatrixPP getInitialGenomeWideMatrixPP(ChromosomeHandler chromosomeHandler) {
        long genomeLength = chromosomeHandler.getChromosomeFromIndex(0).getLength();  // <= whole genome in KB
        int binSize = (int) (genomeLength / 500); // todo
        if (binSize == 0) binSize = 1;
        int nBinsX = (int) (genomeLength / binSize + 1); // todo
        int nBlockColumns = nBinsX / BLOCK_SIZE + 1;
        return new MatrixPP(0, 0, binSize, nBlockColumns, chromosomeHandler, countThreshold, v9DepthBase);
    }

    /**
     * @param file List of files to read
     * @return Matrix with counts in each bin
     */
    private MatrixPP computeWholeGenomeMatrix(String file) throws IOException {

        MatrixPP matrix = getInitialGenomeWideMatrixPP(chromosomeHandler);

        PairIterator iter = null;

        // int belowMapq = 0;
        // int intraFrag = 0;
        // int totalRead = 0;
        // int contig = 0;
        // int hicContact = 0;

        // Create an index the first time through
        try {
            iter = PairIterator.getIterator(file, chromosomeIndexes, chromosomeHandler);
            //ContactFilter filter = getContactFilter();

            while (iter.hasNext()) {
                // totalRead++;
                AlignmentPair pair = iter.next();
                if (pair.isNotContigPair()) {
                    int bp1 = pair.getPos1();
                    int bp2 = pair.getPos2();
                    int chr1 = pair.getChr1();
                    int chr2 = pair.getChr2();

                    int pos1, pos2;
                    if (shouldSkipContact(pair)) continue;
                    pos1 = ContactCleaner.getWholeGenomePosition(chr1, bp1, chromosomeHandler);
                    pos2 = ContactCleaner.getWholeGenomePosition(chr2, bp2, chromosomeHandler);
                    matrix.incrementCount(pos1, pos2, pos1, pos2, pair.getScore(), expectedValueCalculations, tmpDir);
                    // hicContact++;
                }
            }
        } finally {
            if (iter != null) iter.close();
        }
        matrix.parsingComplete();
        return matrix;
    }



    protected void writeBody(String inputFile, Map<Integer, List<Chunk>> mndIndex) throws IOException {

        MatrixPP wholeGenomeMatrix = computeWholeGenomeMatrix(inputFile);
        writeMatrix(wholeGenomeMatrix, losArray, compressor, matrixPositions, -1, false);

        PairIterator iter = PairIterator.getIterator(inputFile, chromosomeIndexes, chromosomeHandler);

        Set<String> writtenMatrices = Collections.synchronizedSet(new HashSet<>());

        int currentChr1 = -1;
        int currentChr2 = -1;
        MatrixPP currentMatrix = null;
        String currentMatrixKey = null;
        ContactCleaner cleaner = new ContactCleaner(chromosomeHandler);

        while (iter.hasNext()) {
            AlignmentPair pair = iter.next();
            // skip pairs that mapped to contigs
            if (pair.isNotContigPair()) {
                if (shouldSkipContact(pair)) continue;
                // Flip pair if needed so chr1 < chr2
                cleaner.updateLatestContact(pair);

                // Randomize position within fragment site

                // only increment if not intraFragment and passes the mapq threshold
                if (cleaner.doesntMatchCurrentBlock(currentChr1, currentChr2)) {
                    // Starting a new matrix
                    if (currentMatrix != null) {
                        currentMatrix.parsingComplete();
                        writeMatrix(currentMatrix, losArray, compressor, matrixPositions, -1, false);
                        writtenMatrices.add(currentMatrixKey);
                        currentMatrix = null;
                        System.gc();
                        //System.out.println("Available memory: " + RuntimeUtils.getAvailableMemory());
                    }

                    // Start the next matrix
                    currentChr1 = cleaner.getChr1();
                    currentChr2 = cleaner.getChr2();
                    currentMatrixKey = currentChr1 + "_" + currentChr2;

                    if (writtenMatrices.contains(currentMatrixKey)) {
                        System.err.println("Error: the chromosome combination " + currentMatrixKey + " appears in multiple blocks");
                        if (outputFile != null) outputFile.deleteOnExit();
                        System.exit(58);
                    }
                    currentMatrix = new MatrixPP(currentChr1, currentChr2, chromosomeHandler, bpBinSizes,
                            countThreshold, v9DepthBase, BLOCK_CAPACITY);
                }
                cleaner.incrementCount(currentMatrix, expectedValueCalculations, tmpDir);
            }
        }

        if (currentMatrix != null) {
            currentMatrix.parsingComplete();
            writeMatrix(currentMatrix, losArray, compressor, matrixPositions, -1, false);
        }

        iter.close();

        masterIndexPosition = losArray[0].getWrittenCount();
    }

    protected boolean shouldSkipContact(AlignmentPair pair) { // static
        int chr1 = pair.getChr1();
        int chr2 = pair.getChr2();
        if (diagonalsOnly && chr1 != chr2) return true;
        if (includedChromosomes != null && chr1 != 0) {
            String c1Name = chromosomeHandler.getChromosomeFromIndex(chr1).getName();
            String c2Name = chromosomeHandler.getChromosomeFromIndex(chr2).getName();
            if (!includedChromosomes.contains(c1Name) || !includedChromosomes.contains(c2Name)) {
                return true;
            }
        }
        if (filter != null) {
            if (!filter.pairTypesAreEqual(pair)) {
                return true;
            }
        }
        int mapq = Math.min(pair.getMapq1(), pair.getMapq2());
        return mapq < mapqThreshold;
    }

    protected void updateMasterIndex(String headerFile) throws IOException {
        try (RandomAccessFile raf = new RandomAccessFile(headerFile, "rw")) {
            // Master index
            raf.getChannel().position(masterIndexPositionPosition);
            BufferedByteWriter buffer = new BufferedByteWriter();
            buffer.putLong(masterIndexPosition);
            raf.write(buffer.getBytes());
            System.out.println("masterIndexPosition: " + masterIndexPosition);
        }
    }

    protected void writeFooter(LittleEndianOutputStream[] los) throws IOException {

        // Index
        BigListOfByteWriters bufferList = new BigListOfByteWriters();

        bufferList.putInt(matrixPositions.size());
        for (Map.Entry<String, IndexEntry> entry : matrixPositions.entrySet()) {
            bufferList.expandBufferIfNeeded(1000);
            bufferList.putNullTerminatedString(entry.getKey());
            bufferList.putLong(entry.getValue().position);
            bufferList.putInt(entry.getValue().size);
        }

        // Vectors - Expected values

        bufferList.expandBufferIfNeeded(1000);
        bufferList.putInt(expectedValueCalculations.size());
        for (Map.Entry<String, ExpectedValueCalculation> entry : expectedValueCalculations.entrySet()) {
            ExpectedValueCalculation ev = entry.getValue();

            ev.computeDensity();

            int binSize = ev.getGridSize();
            HiCZoom.HiCUnit unit = HiCZoom.HiCUnit.BP;

            bufferList.putNullTerminatedString(unit.toString());
            bufferList.putInt(binSize);

            // The density values
            ListOfDoubleArrays expectedValues = ev.getDensityAvg();
            // todo @Suhas to handle buffer overflow
            bufferList.putLong(expectedValues.getLength());
            for (double[] expectedArray : expectedValues.getValues()) {
                bufferList.expandBuffer();
                for (double value : expectedArray) {
                    bufferList.expandBufferIfNeeded(1000000);
                    bufferList.putFloat((float) value);
                }
            }

            // Map of chromosome index -> normalization factor
            Map<Integer, Double> normalizationFactors = ev.getChrScaleFactors();
            bufferList.expandBufferIfNeeded(1000000);
            bufferList.putInt(normalizationFactors.size());
            for (Map.Entry<Integer, Double> normFactor : normalizationFactors.entrySet()) {
                bufferList.putInt(normFactor.getKey());
                bufferList.putFloat(normFactor.getValue().floatValue());
                //System.out.println(normFactor.getKey() + "  " + normFactor.getValue());
            }
        }


        long nBytesV5 = bufferList.getTotalBytes();
        System.out.println("nBytesV5: " + nBytesV5);

        los[0].writeLong(nBytesV5);
        bufferList.writeToOutput(los[0]);
    }

    protected Pair<Map<Long, List<IndexEntry>>, Long> writeMatrix(MatrixPP matrix, LittleEndianOutputStream[] losArray,
                                                                  Deflater compressor, Map<String, IndexEntry> matrixPositions, int chromosomePairIndex, boolean doMultiThreadedBehavior) throws IOException {

        LittleEndianOutputStream los = losArray[0];
        long position = los.getWrittenCount();

        los.writeInt(matrix.getChr1Idx());
        los.writeInt(matrix.getChr2Idx());
        int numResolutions = 0;

        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null) {
                numResolutions++;
            }
        }
        los.writeInt(numResolutions);

        //fos.writeInt(matrix.getZoomData().length);
        for ( int i = 0; i < matrix.getZoomData().length; i++) {
            MatrixZoomDataPP zd = matrix.getZoomData()[i];
            if (zd != null) {
                WriterUtils.writeZoomHeader(zd, los);
            }
        }

        long size = los.getWrittenCount() - position;
        if (chromosomePairIndex > -1) {
            matrixPositions.put("" + chromosomePairIndex, new IndexEntry(position, (int) size));
        } else {
            matrixPositions.put(matrix.getKey(), new IndexEntry(position, (int) size));
        }

        final Map<Long, List<IndexEntry>> localBlockIndexes = new ConcurrentHashMap<>();

        for (int i = 0; i < matrix.getZoomData().length; i++) {
            MatrixZoomDataPP zd = matrix.getZoomData()[i];
            if (zd != null) {
                List<IndexEntry> blockIndex;
                if (doMultiThreadedBehavior) {
                    if (losArray.length > 1) {
                        blockIndex = zd.mergeAndWriteBlocks(losArray, i, matrix.getZoomData().length);
                    } else {
                        blockIndex = zd.mergeAndWriteBlocks(losArray[0], compressor);
                    }
                    localBlockIndexes.put(zd.blockIndexPosition, blockIndex);
                } else {
                    blockIndex = zd.mergeAndWriteBlocks(losArray[0], compressor);
                    updateIndexPositions(blockIndex, losArray, true, outputFile, 0, zd.blockIndexPosition);
                }
            }
        }

        System.out.print(".");
        return new Pair<>(localBlockIndexes, position);
    }

    protected void updateIndexPositions(List<IndexEntry> blockIndex, LittleEndianOutputStream[] losArray, boolean doRestore,
                                        File outputFile, long currentPosition, long blockIndexPosition) throws IOException {

        // Temporarily close output stream.  Remember position
        long losPos = 0;
        if (doRestore) {
            losPos = losArray[0].getWrittenCount();
            losArray[0].close();
        }

        try (RandomAccessFile raf = new RandomAccessFile(outputFile, "rw")) {

            // Block indices
            raf.getChannel().position(blockIndexPosition);

            // Write as little endian
            BufferedByteWriter buffer = new BufferedByteWriter();
            for (IndexEntry aBlockIndex : blockIndex) {
                buffer.putInt(aBlockIndex.id);
                buffer.putLong(aBlockIndex.position + currentPosition);
                buffer.putInt(aBlockIndex.size);
            }
            raf.write(buffer.getBytes());

        }
        if (doRestore) {
            FileOutputStream fos = new FileOutputStream(outputFile, true);
            fos.getChannel().position(losPos);
            losArray[0] = new LittleEndianOutputStream(new BufferedOutputStream(fos, HiCGlobals.bufferSize));
            losArray[0].setWrittenCount(losPos);
        }
    }
}
