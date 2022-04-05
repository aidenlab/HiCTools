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

import hic.HiCGlobals;
import hic.tools.utils.iterators.mnd.ReadPairFilter;
import hic.tools.utils.largelists.BigListOfByteWriters;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.tools.UNIXTools;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.Pair;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.Deflater;

abstract public class HiCFileBuilder {

    protected long masterIndexPositionPosition;
    protected long normVectorIndexPosition;
    protected long normVectorLengthPosition;
    protected Map<String, ExpectedValueCalculation> expectedValueCalculations = Collections.synchronizedMap(new LinkedHashMap<>());

    protected static final int VERSION = 9;
    protected static final int BLOCK_SIZE = 1000;
    public static int BLOCK_CAPACITY = 1000;
    protected final ChromosomeHandler chromosomeHandler;
    protected final File outputFile;
    protected final Map<String, IndexEntry> matrixPositions = new LinkedHashMap<>();
    protected final Deflater compressor = WriterUtils.getDefaultCompressor();
    protected int v9DepthBase = 2;
    protected Map<String, Integer> chromosomeIndexes = new Hashtable<>();
    protected String genomeId;
    protected final LittleEndianOutputStream[] losArray = new LittleEndianOutputStream[1];
    protected long masterIndexPosition;
    protected int countThreshold = 0;
    protected int mapqThreshold = 0;
    protected boolean diagonalsOnly = false;
    protected String statsFileName = null;
    protected String graphFileName = null;
    protected Set<String> includedChromosomes;
    protected ReadPairFilter filter = null;
    protected int[] bpBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000};
    protected int numResolutions;
    protected double hicFileScalingFactor = 1;
    protected File tmpDir = UNIXTools.makeDir(new File("temp_folder"));

    public HiCFileBuilder(File outputFile, String genomeId, ChromosomeHandler chromosomeHandler, double hicFileScalingFactor) {
        this.genomeId = genomeId;
        this.outputFile = outputFile;
        this.chromosomeHandler = chromosomeHandler;
        for (int i = 0; i < chromosomeHandler.size(); i++) {
            chromosomeIndexes.put(chromosomeHandler.getChromosomeFromIndex(i).getName(), i);
        }
        if (hicFileScalingFactor > 0) {
            this.hicFileScalingFactor = hicFileScalingFactor;
        }
        initializeExpectedVectorCalculations();
    }

    protected static void closeLosArray(LittleEndianOutputStream[] los) throws IOException {
        if (los != null && los[0] != null) {
            los[0].close();
        }
    }

    public void setCountThreshold(int countThreshold) {
        this.countThreshold = countThreshold;
    }

    public void setV9DepthBase(int v9DepthBase) {
        if (v9DepthBase > 1 || v9DepthBase < 0) {
            this.v9DepthBase = v9DepthBase;
        }
    }

    public void setMapqThreshold(int mapqThreshold) {
        this.mapqThreshold = mapqThreshold;
    }

    public void setDiagonalsOnly(boolean diagonalsOnly) {
        this.diagonalsOnly = diagonalsOnly;
    }

    public void setIncludedChromosomes(Set<String> includedChromosomes) {
        if (includedChromosomes != null && includedChromosomes.size() > 0) {
            this.includedChromosomes = Collections.synchronizedSet(new HashSet<>());
            for (String name : includedChromosomes) {
                this.includedChromosomes.add(chromosomeHandler.cleanUpName(name));
            }
        }
    }

    public void setGraphFile(String graphFileName) {
        this.graphFileName = graphFileName;
    }

    public void setGenome(String genome) {
        if (genome != null) {
            this.genomeId = genome;
        }
    }

    protected static void updateIndexPositions(List<IndexEntry> blockIndex, LittleEndianOutputStream[] losArray, boolean doRestore,
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

    public void setFilter(ReadPairFilter.Type type) {
        if (type != null) {
            this.filter = new ReadPairFilter(type);
        }
    }


    public void setTmpdir(String tmpDirName) {
        if (tmpDirName != null) {
            this.tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                System.err.println("Tmp directory does not exist: " + tmpDirName);
                if (outputFile != null) outputFile.deleteOnExit();
                System.exit(59);
            }
        }
    }

    public void setStatisticsFile(String statsOption) {
        statsFileName = statsOption;
    }

    private void initializeExpectedVectorCalculations() {
        expectedValueCalculations.clear();
        for (int bBinSize : bpBinSizes) {
            ExpectedValueCalculation calc = new ExpectedValueCalculation(chromosomeHandler, bBinSize, NormalizationHandler.NONE);
            String key = "BP_" + bBinSize;
            expectedValueCalculations.put(key, calc);
        }
    }

    public void setResolutions(List<String> resolutions) {
        if (resolutions != null) {
            ArrayList<Integer> bpResolutions = new ArrayList<>();
            for (String str : resolutions) {
                try {
                    int myInt = Integer.parseInt(str);
                    bpResolutions.add(myInt);
                } catch (NumberFormatException exception) {
                    System.err.println("Resolution improperly formatted. It must be in the form of a number, such as 1000000 for 1 million BP");
                    System.err.println("Fragment resolutions are deprecated and require using an older jar.");
                    System.exit(1);
                }
            }

            boolean resolutionsSet = false;
            if (bpResolutions.size() > 0) {
                resolutionsSet = true;
                Collections.sort(bpResolutions);
                Collections.reverse(bpResolutions);
                int[] bps = new int[bpResolutions.size()];
                for (int i = 0; i < bps.length; i++) {
                    bps[i] = bpResolutions.get(i);
                }
                bpBinSizes = bps;
                // reset expected vectors based on new resolutions
                initializeExpectedVectorCalculations();
            } else {
                bpBinSizes = new int[0];
            }
            if (!resolutionsSet) {
                System.err.println("No valid resolutions sent in");
                System.exit(1);
            }
        }
    }

    protected LittleEndianOutputStream[] initializeLosArrays(String headerFile, String footerFile) {
        try {
            losArray[0] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(headerFile), HiCGlobals.bufferSize));
            if (footerFile.equalsIgnoreCase(headerFile)) {
                return losArray;
            } else {
                LittleEndianOutputStream[] losFooter = new LittleEndianOutputStream[1];
                losFooter[0] = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(footerFile), HiCGlobals.bufferSize));
                return losFooter;
            }
        } catch (Exception e) {
            System.err.println("Unable to write to " + outputFile);
            System.exit(70);
        }
        return null;
    }

    protected StringBuilder readFileIntoString(String fileName, String description) {
        if (fileName != null) {
            try (FileInputStream is = new FileInputStream(fileName)) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(is), HiCGlobals.bufferSize);
                StringBuilder sb = new StringBuilder();
                String nextLine;
                while ((nextLine = reader.readLine()) != null) {
                    sb.append(nextLine).append("\n");
                }
                return sb;
            } catch (IOException e) {
                System.err.println("Error while reading " + description + " file: " + e);
                return null;
            }
        }
        return null;
    }

    protected MatrixPP getInitialGenomeWideMatrixPP(ChromosomeHandler chromosomeHandler) {
        long genomeLength = chromosomeHandler.getChromosomeFromIndex(0).getLength();  // <= whole genome in KB
        int binSize = (int) (genomeLength / 500); // todo
        if (binSize == 0) binSize = 1;
        int nBinsX = (int) (genomeLength / binSize + 1); // todo
        int nBlockColumns = nBinsX / BLOCK_SIZE + 1;
        return new MatrixPP(0, 0, binSize, nBlockColumns, chromosomeHandler, countThreshold, v9DepthBase);
    }

    protected void writeHeader() throws IOException {
        System.out.println("Start preprocess");
        System.out.println("Writing header");
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

        StringBuilder stats = readFileIntoString(statsFileName, "stats");
        StringBuilder graphs = readFileIntoString(graphFileName, "graphs");
        StringBuilder hicFileScaling = new StringBuilder().append(hicFileScalingFactor);

        // Attribute dictionary
        int nAttributes = 2;
        if (stats != null) nAttributes += 1;
        if (graphs != null) nAttributes += 1;
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
        { // hicFileScaling always written
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

    protected void writeFooter(LittleEndianOutputStream[] los) throws IOException {
        System.out.println();
        System.out.println("Writing footer");
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


        long nBytesV5 = bufferList.getBytesWritten();
        System.out.println("nBytesV5: " + nBytesV5);

        los[0].writeLong(nBytesV5);
        bufferList.writeToOutput(los[0]);
    }

    protected Pair<Map<Long, List<IndexEntry>>, Long> writeMatrix(MatrixPP matrix, LittleEndianOutputStream[] losArray,
                                                                  Deflater compressor, Map<String, IndexEntry> matrixPositions,
                                                                  int chromosomePairIndex, boolean doMultiThreadedBehavior) throws IOException {

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
        for (int i = 0; i < matrix.getZoomData().length; i++) {
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
}
