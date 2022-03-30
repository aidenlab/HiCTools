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
import hic.tools.clt.CommandLineParser.Alignment;
import hic.tools.utils.cleaner.ContactCleaner;
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
public class Preprocessor {
    protected static final int VERSION = 9;
    protected static final int BLOCK_SIZE = 1000;
    protected int v9DepthBase = 2;
    protected final ChromosomeHandler chromosomeHandler;
    protected Map<String, Integer> chromosomeIndexes;
    protected final File outputFile;
    protected final Map<String, IndexEntry> matrixPositions;
    protected String genomeId;
    protected final Deflater compressor;
    protected LittleEndianOutputStream[] losArray = new LittleEndianOutputStream[1];
    protected long masterIndexPosition;
    protected int countThreshold = 0;
    protected int mapqThreshold = 0;
    protected boolean diagonalsOnly = false;
    protected String statsFileName = null;
    protected String graphFileName = null;
    protected String expectedVectorFile = null;
    protected Set<String> randomizeFragMapFiles = null;
    protected Set<String> includedChromosomes;
    protected Alignment alignmentFilter;
    public static int BLOCK_CAPACITY = 1000;
    
    // Base-pair resolutions
    protected int[] bpBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000};
    
    // Fragment resolutions
    protected int[] fragBinSizes = {500, 200, 100, 50, 20, 5, 2, 1};

    // number of resolutions
    protected int numResolutions = bpBinSizes.length + fragBinSizes.length;

    // hic scaling factor value
    protected double hicFileScalingFactor = 1;
    
    protected Long normVectorIndex = 0L, normVectorLength = 0L;
    
    /**
     * The position of the field containing the masterIndex position
     */
    protected long masterIndexPositionPosition;
    protected long normVectorIndexPosition;
    protected long normVectorLengthPosition;
    protected Map<String, ExpectedValueCalculation> expectedValueCalculations;
    protected File tmpDir;
    
    public Preprocessor(File outputFile, String genomeId, ChromosomeHandler chromosomeHandler, double hicFileScalingFactor) {
        this.genomeId = genomeId;
        this.outputFile = outputFile;
        this.matrixPositions = new LinkedHashMap<>();

        this.chromosomeHandler = chromosomeHandler;
        chromosomeIndexes = new Hashtable<>();
        for (int i = 0; i < chromosomeHandler.size(); i++) {
            chromosomeIndexes.put(chromosomeHandler.getChromosomeFromIndex(i).getName(), i);
        }

        compressor = WriterUtils.getDefaultCompressor();

        this.tmpDir = null;  // TODO -- specify this

        if (hicFileScalingFactor > 0) {
            this.hicFileScalingFactor = hicFileScalingFactor;
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

    public void setExpectedVectorFile(String expectedVectorFile) {
        this.expectedVectorFile = expectedVectorFile;
    }

    public void setGraphFile(String graphFileName) {
        this.graphFileName = graphFileName;
    }

    public void setGenome(String genome) {
        if (genome != null) {
            this.genomeId = genome;
        }
    }

    public void setResolutions(List<String> resolutions) {
        if (resolutions != null) {
            ArrayList<Integer> fragResolutions = new ArrayList<>();
            ArrayList<Integer> bpResolutions = new ArrayList<>();

            for (String str : resolutions) {
                boolean fragment = false;
                int index = str.indexOf("f");
                if (index != -1) {
                    str = str.substring(0, index);
                    fragment = true;
                }
                Integer myInt = null;
                try {
                    myInt = Integer.valueOf(str);
                } catch (NumberFormatException exception) {
                    System.err.println("Resolution improperly formatted.  It must be in the form of a number, such as 1000000 for 1M bp,");
                    System.err.println("or a number followed by 'f', such as 25f for 25 fragment");
                    System.exit(1);
                }
                if (fragment) fragResolutions.add(myInt);
                else          bpResolutions.add(myInt);
            }

            boolean resolutionsSet = false;
            if (fragResolutions.size() > 0) {
                resolutionsSet = true;
                Collections.sort(fragResolutions);
                Collections.reverse(fragResolutions);
                int[] frags = new int[fragResolutions.size()];
                for (int i=0; i<frags.length; i++){
                    frags[i] = fragResolutions.get(i);
                }
                fragBinSizes = frags;
            }
            else {
                fragBinSizes = new int[0];
            }
            if (bpResolutions.size() > 0) {
                resolutionsSet = true;
                Collections.sort(bpResolutions);
                Collections.reverse(bpResolutions);
                int[] bps = new int[bpResolutions.size()];
                for (int i = 0; i < bps.length; i++) {
                    bps[i] = bpResolutions.get(i);
                }
                bpBinSizes = bps;
            }
            else {
                bpBinSizes = new int[0];
            }
            if (!resolutionsSet) {
                System.err.println("No valid resolutions sent in");
                System.exit(1);
            }
        }
    }

    public void setAlignmentFilter(Alignment al) {
        this.alignmentFilter = al;
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

            if (expectedVectorFile == null) {
                expectedValueCalculations = Collections.synchronizedMap(new LinkedHashMap<>());
                for (int bBinSize : bpBinSizes) {
                    ExpectedValueCalculation calc = new ExpectedValueCalculation(chromosomeHandler, bBinSize, NormalizationHandler.NONE);
                    String key = "BP_" + bBinSize;
                    expectedValueCalculations.put(key, calc);
                }
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
        int nFragRes = 0;
        los.writeInt(nFragRes);

        numResolutions = nBpRes + nFragRes;
    }

    protected MatrixPP getInitialGenomeWideMatrixPP(ChromosomeHandler chromosomeHandler) {
        long genomeLength = chromosomeHandler.getChromosomeFromIndex(0).getLength();  // <= whole genome in KB
        int binSize = (int) (genomeLength / 500); // todo
        if (binSize == 0) binSize = 1;
        int nBinsX = (int) (genomeLength / binSize + 1); // todo
        int nBlockColumns = nBinsX / BLOCK_SIZE + 1;
        return new MatrixPP(0, 0, binSize, nBlockColumns, chromosomeHandler, countThreshold, v9DepthBase);
    }

    protected boolean alignmentsAreEqual(Alignment alignment, Alignment alignmentStandard) {
        if (alignment == alignmentStandard) {
            return true;
        }
        if (alignmentStandard == Alignment.TANDEM) {
            return alignment == Alignment.LL || alignment == Alignment.RR;
        }

        return false;
    }

    /**
     * @param file List of files to read
     * @return Matrix with counts in each bin
     * @throws IOException
     */
    private MatrixPP computeWholeGenomeMatrix(String file) throws IOException {

        MatrixPP matrix = getInitialGenomeWideMatrixPP(chromosomeHandler);

        PairIterator iter = null;

        //int belowMapq = 0;
        //int intraFrag = 0;
        int totalRead = 0;
        int contig = 0;
        int hicContact = 0;

        // Create an index the first time through
        try {
            iter = PairIterator.getIterator(file, chromosomeIndexes, chromosomeHandler);
            //ContactFilter filter = getContactFilter();

            while (iter.hasNext()) {
                totalRead++;
                AlignmentPair pair = iter.next();
                if (pair.isContigPair()) {
                    contig++;
                } else {
                    int bp1 = pair.getPos1();
                    int bp2 = pair.getPos2();
                    int chr1 = pair.getChr1();
                    int chr2 = pair.getChr2();

                    int pos1, pos2;
                    if (shouldSkipContact(pair)) continue;
                    pos1 = ContactCleaner.getWholeGenomePosition(chr1, bp1, chromosomeHandler);
                    pos2 = ContactCleaner.getWholeGenomePosition(chr2, bp2, chromosomeHandler);
                    matrix.incrementCount(pos1, pos2, pos1, pos2, pair.getScore(), expectedValueCalculations, tmpDir);
                    hicContact++;
                }
            }
        } finally {
            if (iter != null) iter.close();
        }
        matrix.parsingComplete();
        return matrix;
    }

    protected static Alignment calculateAlignment(AlignmentPair pair) {

        if (pair.getStrand1() == pair.getStrand2()) {
            if (pair.getStrand1()) {
                return Alignment.RR;
            } else {
                return Alignment.LL;
            }
        } else if (pair.getStrand1()) {
            if (pair.getPos1() < pair.getPos2()) {
                return Alignment.INNER;
            } else {
                return Alignment.OUTER;
            }
        } else {
            if (pair.getPos1() < pair.getPos2()) {
                return Alignment.OUTER;
            } else {
                return Alignment.INNER;
            }
        }
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
            if (!pair.isContigPair()) {
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
        if (alignmentFilter != null && !alignmentsAreEqual(calculateAlignment(pair), alignmentFilter)) {
            return true;
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
        List<BufferedByteWriter> bufferList = new ArrayList<>();
        bufferList.add(new BufferedByteWriter());
        bufferList.get(0).putInt(matrixPositions.size());
        for (Map.Entry<String, IndexEntry> entry : matrixPositions.entrySet()) {
            if (Integer.MAX_VALUE - bufferList.get(bufferList.size()-1).bytesWritten() < 1000) {
                bufferList.add(new BufferedByteWriter());
            }
            bufferList.get(bufferList.size()-1).putNullTerminatedString(entry.getKey());
            bufferList.get(bufferList.size()-1).putLong(entry.getValue().position);
            bufferList.get(bufferList.size()-1).putInt(entry.getValue().size);
        }

        // Vectors  (Expected values,  other).
        /***  NEVA ***/
        if (expectedVectorFile == null) {
            if (Integer.MAX_VALUE - bufferList.get(bufferList.size()-1).bytesWritten() < 1000) {
                bufferList.add(new BufferedByteWriter());
            }
            bufferList.get(bufferList.size()-1).putInt(expectedValueCalculations.size());
            for (Map.Entry<String, ExpectedValueCalculation> entry : expectedValueCalculations.entrySet()) {
                ExpectedValueCalculation ev = entry.getValue();
    
                ev.computeDensity();

                int binSize = ev.getGridSize();
                HiCZoom.HiCUnit unit = HiCZoom.HiCUnit.BP;

                bufferList.get(bufferList.size()-1).putNullTerminatedString(unit.toString());
                bufferList.get(bufferList.size()-1).putInt(binSize);
    
                // The density values
                ListOfDoubleArrays expectedValues = ev.getDensityAvg();
                // todo @Suhas to handle buffer overflow
                bufferList.get(bufferList.size()-1).putLong(expectedValues.getLength());
                for (double[] expectedArray : expectedValues.getValues()) {
                    bufferList.add(new BufferedByteWriter());
                    for (double value : expectedArray) {
                        if (Integer.MAX_VALUE - bufferList.get(bufferList.size()-1).bytesWritten() < 1000000) {
                            bufferList.add(new BufferedByteWriter());
                        }
                        bufferList.get(bufferList.size()-1).putFloat( (float) value);
                    }
                }
    
                // Map of chromosome index -> normalization factor
                Map<Integer, Double> normalizationFactors = ev.getChrScaleFactors();
                if (Integer.MAX_VALUE - bufferList.get(bufferList.size()-1).bytesWritten() < 1000000) {
                    bufferList.add(new BufferedByteWriter());
                }
                bufferList.get(bufferList.size()-1).putInt(normalizationFactors.size());
                for (Map.Entry<Integer, Double> normFactor : normalizationFactors.entrySet()) {
                    bufferList.get(bufferList.size()-1).putInt(normFactor.getKey());
                    bufferList.get(bufferList.size()-1).putFloat(normFactor.getValue().floatValue());
                    //System.out.println(normFactor.getKey() + "  " + normFactor.getValue());
                }
            }
        }
        else {
            // read in expected vector file. to get # of resolutions, might have to read twice.

            int count=0;
            try (Reader reader = new FileReader(expectedVectorFile);
                 BufferedReader bufferedReader = new BufferedReader(reader)) {

                String line;
                while ((line = bufferedReader.readLine()) != null) {
                    if (line.startsWith("fixedStep"))
                        count++;
                    if (line.startsWith("variableStep")) {
                        System.err.println("Expected vector file must be in wiggle fixedStep format");
                        System.exit(19);
                    }
                }
            }
            bufferList.get(bufferList.size()-1).putInt(count);
            try (Reader reader = new FileReader(expectedVectorFile);
                 BufferedReader bufferedReader = new BufferedReader(reader)) {

                String line;
                while ((line = bufferedReader.readLine()) != null) {
                    if (line.startsWith("fixedStep")) {
                        String[] words = line.split("\\s+");
                        for (String str:words){
                            if (str.contains("chrom")){
                                String[] chrs = str.split("=");

                            }
                        }
                    }
                }
            }
        }
        long nBytesV5 = 0;
        for (int i = 0; i<bufferList.size(); i++) {
            nBytesV5 += bufferList.get(i).getBytes().length;
        }
        System.out.println("nBytesV5: " + nBytesV5);

        los[0].writeLong(nBytesV5);
        for (int i = 0; i<bufferList.size(); i++) {
            los[0].write(bufferList.get(i).getBytes());
        }
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
                List<IndexEntry> blockIndex = null;
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

        RandomAccessFile raf = null;
        try {
            raf = new RandomAccessFile(outputFile, "rw");

            // Block indices
            long pos = blockIndexPosition;
            raf.getChannel().position(pos);

            // Write as little endian
            BufferedByteWriter buffer = new BufferedByteWriter();
            for (IndexEntry aBlockIndex : blockIndex) {
                buffer.putInt(aBlockIndex.id);
                buffer.putLong(aBlockIndex.position + currentPosition);
                buffer.putInt(aBlockIndex.size);
            }
            raf.write(buffer.getBytes());

        } finally {
            if (raf != null) raf.close();
        }
        if (doRestore) {
            FileOutputStream fos = new FileOutputStream(outputFile, true);
            fos.getChannel().position(losPos);
            losArray[0] = new LittleEndianOutputStream(new BufferedOutputStream(fos, HiCGlobals.bufferSize));
            losArray[0].setWrittenCount(losPos);
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
}
