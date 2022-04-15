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

import htsjdk.tribble.util.LittleEndianInputStream;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.depth.V9Depth;
import javastraw.tools.ParallelizationTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;

import java.awt.*;
import java.io.*;
import java.util.List;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.Deflater;

public class MatrixZoomDataPP {
    final Set<Integer> blockNumbers;  // The only reason for this is to get a count
    final ConcurrentHashMap<Integer, Integer> blockNumRecords;
    final List<File> tmpFiles;
    final Map<Integer, Map<File, Long>> tmpFilesByBlockNumber;
    private final Chromosome chr1;  // Redundant, but convenient    BinDatasetReader
    private final Chromosome chr2;  // Redundant, but convenient
    private final int zoom;
    private final int binSize;              // bin size in bp
    private final int blockBinCount;        // block size in bins
    private final int blockColumnCount;     // number of block columns
    private final LinkedHashMap<Integer, BlockPP> blocks;
    private final int countThreshold;
    long blockIndexPosition;
    private double sum = 0;
    private double numRecords = 0;
    private double cellCount = 0;
    private double percent5;
    private double percent95;
    private int BLOCK_CAPACITY = 1000;
    private final V9Depth v9Depth;

    /**
     * Representation of MatrixZoomData used for preprocessing
     *
     * @param chr1             index of first chromosome  (x-axis)
     * @param chr2             index of second chromosome
     * @param binSize          size of each grid bin in bp
     * @param blockColumnCount number of block columns
     * @param zoom             integer zoom (resolution) level index.  TODO Is this needed?
     */
    MatrixZoomDataPP(Chromosome chr1, Chromosome chr2, int binSize, int blockColumnCount, int zoom, int countThreshold, int v9BaseDepth) {
        this.tmpFiles = new ArrayList<>();
        this.tmpFilesByBlockNumber = new ConcurrentHashMap<>();
        this.blockNumbers = Collections.synchronizedSet(new HashSet<>(BLOCK_CAPACITY));
        this.blockNumRecords = new ConcurrentHashMap<>(BLOCK_CAPACITY);
        this.countThreshold = countThreshold;

        this.chr1 = chr1;
        this.chr2 = chr2;
        this.binSize = binSize;
        this.blockColumnCount = blockColumnCount;
        this.zoom = zoom;

        // Get length in proper units
        Chromosome longChr = chr1.getLength() > chr2.getLength() ? chr1 : chr2;
        long len = longChr.getLength();

        int nBinsX = (int) (len / binSize + 1);

        blockBinCount = nBinsX / blockColumnCount + 1;
        blocks = new LinkedHashMap<>(blockColumnCount);
        v9Depth = V9Depth.setDepthMethod(v9BaseDepth, blockBinCount);
    }

    MatrixZoomDataPP(Chromosome chr1, Chromosome chr2, int binSize, int blockColumnCount, int zoom,
                     int countThreshold, int v9BaseDepth, int BLOCK_CAPACITY) {
        this.tmpFiles = new ArrayList<>();
        this.tmpFilesByBlockNumber = new ConcurrentHashMap<>();
        this.BLOCK_CAPACITY = BLOCK_CAPACITY;
        this.blockNumbers = Collections.synchronizedSet(new HashSet<>(BLOCK_CAPACITY));
        this.blockNumRecords = new ConcurrentHashMap<>(BLOCK_CAPACITY);
        this.countThreshold = countThreshold;

        this.chr1 = chr1;
        this.chr2 = chr2;
        this.binSize = binSize;
        this.blockColumnCount = blockColumnCount;
        this.zoom = zoom;

        // Get length in proper units
        Chromosome longChr = chr1.getLength() > chr2.getLength() ? chr1 : chr2;
        long len = longChr.getLength();

        int nBinsX = (int) (len / binSize + 1);

        blockBinCount = nBinsX / blockColumnCount + 1;
        blocks = new LinkedHashMap<>(blockColumnCount);
        v9Depth = V9Depth.setDepthMethod(v9BaseDepth, blockBinCount);
    }

    double getSum() {
        return sum;
    }

    double getOccupiedCellCount() {
        return cellCount;
    }

    double getPercent95() {
        return percent95;
    }

    double getPercent5() {
        return percent5;
    }

    int getBinSize() {
        return binSize;
    }

    int getZoom() {
        return zoom;
    }

    int getBlockBinCount() {
        return blockBinCount;
    }

    int getBlockColumnCount() {
        return blockColumnCount;
    }

    /**
     * Increment the count for the bin represented by the GENOMIC position (pos1, pos2)
     */
    void incrementCount(int pos1, int pos2, float score, Map<String, ExpectedValueCalculation> expectedValueCalculations,
                        File tmpDir) throws IOException {
        if (pos1 < 0 || pos2 < 0) return;
        sum += score;
        int xBin = pos1 / binSize;
        int yBin = pos2 / binSize;
        commonIncrementCount(xBin, yBin, score, expectedValueCalculations, tmpDir);
    }

    private void commonIncrementCount(int xBin0, int yBin0, float score,
                                      Map<String, ExpectedValueCalculation> expectedValueCalculations,
                                      File tmpDir) throws IOException {
        int blockNumber;
        int xBin = xBin0;
        int yBin = yBin0;

        // Intra chromosome -- we'll store lower diagonal only
        if (chr1.equals(chr2)) {
            xBin = Math.min(xBin0, yBin0);
            yBin = Math.max(xBin0, yBin0);

            if (xBin != yBin) {
                sum += score;  // <= count for mirror cell.
            }

            if (expectedValueCalculations != null) {
                String evKey = "BP_" + binSize;
                ExpectedValueCalculation ev = expectedValueCalculations.get(evKey);
                if (ev != null) {
                    ev.addDistance(chr1.getIndex(), xBin, yBin, score);
                }
            }

            //compute intra chromosomal block number (version 9 and up)
            int depth = v9Depth.getDepth(xBin, yBin);
            int positionAlongDiagonal = ((xBin + yBin) / 2 / blockBinCount);
            blockNumber = depth * blockColumnCount + positionAlongDiagonal;
        } else {
            // compute interchromosomal block number (version 9 and up, first block is zero)
            int blockCol = xBin0 / blockBinCount;
            int blockRow = yBin0 / blockBinCount;
            blockNumber = blockColumnCount * blockRow + blockCol;
        }

        BlockPP block = blocks.get(blockNumber);
        if (block == null) {
            block = new BlockPP(blockNumber);
            blocks.put(blockNumber, block);
        }
        block.incrementCount(xBin, yBin, score);

        if (blocks.size() > BLOCK_CAPACITY) {
            File tmpfile;
            if (tmpDir == null) {
                tmpfile = File.createTempFile("blocks", "bin");
            } else {
                tmpfile = File.createTempFile("blocks", "bin", tmpDir);
            }
            dumpBlocks(tmpfile);
            tmpFiles.add(tmpfile);
            tmpfile.deleteOnExit();
        }
    }

    static Map<Point, Float> readContactRecordsToMap(int nRecords, byte[] bytes) throws IOException {
        ByteArrayInputStream bis = new ByteArrayInputStream(bytes);
        LittleEndianInputStream lis = new LittleEndianInputStream(bis);

        Map<Point, Float> contactRecordMap = new HashMap<>(nRecords);
        for (int i = 0; i < nRecords; i++) {
            int x = lis.readInt();
            int y = lis.readInt();
            float v = lis.readFloat();
            contactRecordMap.put(new Point(x, y), v);
        }
        try {
            lis.close();
            bis.close();
        } catch (Exception e) {
            System.err.println("Error cleanup contact record to map = " + e.getLocalizedMessage());
        }
        return contactRecordMap;
    }

    /**
     * todo should this be synchronized?
     *
     * @param data
     * @param compressor
     * @return
     */
    private static byte[] compress(byte[] data, Deflater compressor) {

        // Give the compressor the data to compress
        compressor.reset();
        compressor.setInput(data);
        compressor.finish();

        // Create an expandable byte array to hold the compressed data.
        // You cannot use an array that's the same size as the orginal because
        // there is no guarantee that the compressed data will be smaller than
        // the uncompressed data.
        ByteArrayOutputStream bos = new ByteArrayOutputStream(data.length);

        // Compress the data
        byte[] buf = new byte[1024];
        while (!compressor.finished()) {
            int count = compressor.deflate(buf);
            bos.write(buf, 0, count);
        }
        try {
            bos.close();
        } catch (IOException e) {
            System.err.println("Error closing ByteArrayOutputStream");
            e.printStackTrace();
        }

        return bos.toByteArray();
    }


    // Merge and write out blocks one at a time.
    protected List<IndexEntry> mergeAndWriteBlocks(LittleEndianOutputStream los, Deflater compressor) throws IOException {
        DownsampledDoubleArrayList sampledData = new DownsampledDoubleArrayList(10000, 10000);

        List<BlockQueue> activeList = new ArrayList<>();

        // Initialize queues -- first whatever is left over in memory
        if (blocks.size() > 0) {
            BlockQueue bqInMem = new BlockQueueMem(blocks.values());
            activeList.add(bqInMem);
        }
        // Now from files
        for (File file : tmpFiles) {
            BlockQueue bq = new BlockQueueFB(file);
            if (bq.getBlock() != null) {
                activeList.add(bq);
            }
        }

        if (activeList.size() == 0) {
            throw new RuntimeException("No reads in Hi-C contact matrices. This could be because the MAPQ filter is set too high (-q) or because all reads map to the same fragment.");
        }

        List<IndexEntry> indexEntries = new ArrayList<>();
        do {
            activeList.sort(Comparator.comparingInt(o -> o.getBlock().getNumber()));

            BlockQueue topQueue = activeList.get(0);
            BlockPP currentBlock = topQueue.getBlock();
            topQueue.advance();
            int num = currentBlock.getNumber();


            for (int i = 1; i < activeList.size(); i++) {
                BlockQueue blockQueue = activeList.get(i);
                BlockPP block = blockQueue.getBlock();
                if (block.getNumber() == num) {
                    currentBlock.merge(block);
                    blockQueue.advance();
                }
            }

            activeList.removeIf(blockQueue -> blockQueue.getBlock() == null);

            // Output block
            long position = los.getWrittenCount();
            writeBlock(currentBlock, sampledData, los, compressor);
            long size = los.getWrittenCount() - position;

            indexEntries.add(new IndexEntry(num, position, (int) size));

        } while (activeList.size() > 0);


        for (File f : tmpFiles) {
            boolean result = f.delete();
            if (!result) {
                System.out.println("Error while deleting file");
            }
        }

        computeStats(sampledData);

        return indexEntries;
    }

    // Merge and write out blocks multithreaded.
    protected List<IndexEntry> mergeAndWriteBlocks(LittleEndianOutputStream[] losArray, int whichZoom, int numResolutions) {
        DownsampledDoubleArrayList sampledData = new DownsampledDoubleArrayList(10000, 10000);
        Integer[] sortedBlockNumbers = new Integer[blockNumbers.size()];
        blockNumbers.toArray(sortedBlockNumbers);
        Arrays.sort(sortedBlockNumbers);
        Map<Integer, BlockPP> threadSafeBlocks = new ConcurrentHashMap<>(blocks.size());
        threadSafeBlocks.putAll(blocks);
        int numCPUThreads = (losArray.length - 1) / numResolutions;

        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        Map<Integer, Long> blockChunkSizes = new ConcurrentHashMap<>(numCPUThreads);
        Map<Integer, List<IndexEntry>> chunkBlockIndexes = new ConcurrentHashMap<>(numCPUThreads);

        int startBlock =0, endBlock = 0;
        for (int threadNum = 0; threadNum < numCPUThreads; threadNum++) {
            final int whichLos = numCPUThreads * whichZoom + threadNum;
            final int numOfRecordsPerThread = 2 * (int) Math.floor(numRecords / numCPUThreads);
            final int maxNumOfBlocksPerThread = (int) Math.floor((double) sortedBlockNumbers.length / numCPUThreads);
            if (threadNum > 0) {
                startBlock = endBlock;
            }
            int numOfRecords = 0;
            for (int i = startBlock; i < sortedBlockNumbers.length; i++) {
                numOfRecords += blockNumRecords.get(sortedBlockNumbers[i]);
                if (numOfRecords > numOfRecordsPerThread || i - startBlock > maxNumOfBlocksPerThread) {
                    endBlock = i; // i always less than sortedBlockNumbers.length from for loop
                    //endBlock = Math.min(i, sortedBlockNumbers.length);
                    break;
                }
            }
            if (threadNum + 1 == numCPUThreads && endBlock < sortedBlockNumbers.length) {
                endBlock = sortedBlockNumbers.length;
            }
            //System.err.println(binSize + " " + blockNumbers.size() + " " + sortedBlockNumbers.length + " " + startBlock + " " + endBlock);
            if (startBlock >= endBlock) {
                blockChunkSizes.put(threadNum, (long) 0);
                continue;
            }
            final Integer[] threadBlocks = Arrays.copyOfRange(sortedBlockNumbers, startBlock, endBlock);
            List<IndexEntry> indexEntries = new ArrayList<>();
            Runnable worker = () -> {
                try {
                    writeBlockChunk(threadBlocks, threadSafeBlocks, losArray, whichLos, indexEntries, sampledData);
                } catch (Exception e) {
                    e.printStackTrace();
                }

                chunkBlockIndexes.put(whichLos, indexEntries);
            };
            executor.execute(worker);
        }

        ParallelizationTools.shutDownAndWaitUntilDone(executor, 1000);

        long adjust = 0;
        for (int i = 0; i < losArray.length; i++) {
            blockChunkSizes.put(i, losArray[i].getWrittenCount());
            if (i < numCPUThreads * whichZoom) {
                adjust += blockChunkSizes.get(i);
            }
        }
        List<IndexEntry> finalIndexEntries = new ArrayList<>();
        for (int i = numCPUThreads * whichZoom; i < numCPUThreads * (whichZoom + 1); i++) {
            adjust += blockChunkSizes.get(i);
            if (chunkBlockIndexes.get(i) != null) {
                for (int j = 0; j < chunkBlockIndexes.get(i).size(); j++) {
                    finalIndexEntries.add(new IndexEntry(chunkBlockIndexes.get(i).get(j).id, chunkBlockIndexes.get(i).get(j).position + adjust,
                            chunkBlockIndexes.get(i).get(j).size));
                }
            }

        }

        for (File f : tmpFiles) {
            boolean result = f.delete();
            if (!result) {
                System.out.println("Error while deleting file");
            }
        }

        computeStats(sampledData);
        return finalIndexEntries;
    }

    private void writeBlockChunk(Integer[] threadBlocks, Map<Integer, BlockPP> threadSafeBlocks,LittleEndianOutputStream[] losArray,
                                 int threadNum, List<IndexEntry> indexEntries, DownsampledDoubleArrayList sampledData ) throws IOException{
        Deflater compressor = new Deflater();
        compressor.setLevel(Deflater.DEFAULT_COMPRESSION);
        //System.err.println(threadBlocks.length);
        for (int i = 0; i < threadBlocks.length; i++) {
            BlockPP currentBlock = null;
            int num = threadBlocks[i];
            if (threadSafeBlocks.get(num) != null ){
                currentBlock = threadSafeBlocks.get(num);
                if (tmpFilesByBlockNumber.get(num) != null) {
                    for (Map.Entry<File, Long> entry : tmpFilesByBlockNumber.get(num).entrySet()) {
                        readAndMerge(currentBlock, entry);
                    }
                }
            } else if (tmpFilesByBlockNumber.get(num) != null) {
                Iterator<Map.Entry<File, Long>> iter = tmpFilesByBlockNumber.get(num).entrySet().iterator();
                if (iter.hasNext()) {
                    Map.Entry<File, Long> firstEntry = iter.next();
                    currentBlock = readTmpBlock(firstEntry.getKey(), firstEntry.getValue());
                    if (currentBlock != null) {
                        while (iter.hasNext()) {
                            readAndMerge(currentBlock, iter.next());
                        }
                    }
                }
            }

            if (currentBlock != null) {
                long position = losArray[threadNum + 1].getWrittenCount();
                writeBlock(currentBlock, sampledData, losArray[threadNum + 1], compressor);
                long size = losArray[threadNum + 1].getWrittenCount() - position;
                indexEntries.add(new IndexEntry(num, position, (int) size));
            }
        }
    }

    private void readAndMerge(BlockPP currentBlock, Map.Entry<File, Long> entry) throws IOException {
        BlockPP tmpBlock = readTmpBlock(entry.getKey(), entry.getValue());
        if (tmpBlock != null) {
            currentBlock.merge(tmpBlock);
        }
    }

    /**
     * Dump the blocks calculated so far to a temporary file
     *
     * @param file File to write to
     */
    private void dumpBlocks(File file) throws IOException {
        try (LittleEndianOutputStream los = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(file), 4194304))) {

            List<BlockPP> blockList = new ArrayList<>(blocks.values());
            blockList.sort(Comparator.comparingInt(BlockPP::getNumber));

            for (BlockPP b : blockList) {
                blocks.remove(b.getNumber());
                int number = addToBlockAndRecordsSets(b);

                if (tmpFilesByBlockNumber.get(number) == null) {
                    tmpFilesByBlockNumber.put(number, new ConcurrentHashMap<>());
                }
                tmpFilesByBlockNumber.get(number).put(file, los.getWrittenCount());

                los.writeInt(number);
                Map<Point, Float> records = b.getContactRecordMap();

                los.writeInt(records.size());
                for (Map.Entry<Point, Float> entry : records.entrySet()) {
                    Point point = entry.getKey();
                    Float count = entry.getValue();

                    los.writeInt(point.x);
                    los.writeInt(point.y);
                    los.writeFloat(count);
                }
            }

            blocks.clear();
        }
    }

    private BlockPP readTmpBlock(File file, long filePosition) throws IOException {
        if (filePosition >= file.length()) {
            return null;
        }

        try (FileInputStream fis = new FileInputStream(file)) {
            fis.getChannel().position(filePosition);

            LittleEndianInputStream lis = new LittleEndianInputStream(fis);
            int blockNumber = lis.readInt();
            int nRecords = lis.readInt();

            byte[] bytes = new byte[nRecords * 12];
            int len = bytes.length;
            int n = 0;
            while (n < len) {
                int count = fis.read(bytes, n, len - n);
                if (count < 0)
                    throw new EOFException();
                n += count;
            }

            return new BlockPP(blockNumber, readContactRecordsToMap(nRecords, bytes));
        }
    }

    private void computeStats(DownsampledDoubleArrayList sampledData) {
        DescriptiveStatistics stats = new DescriptiveStatistics(sampledData.toArray());
        this.percent5 = stats.getPercentile(5);
        this.percent95 = stats.getPercentile(95);
    }

    void parsingComplete() {
        for (BlockPP block : blocks.values()) {
            addToBlockAndRecordsSets(block);
        }
    }

    /**
     * used by multithreaded code
     * @param otherMatrixZoom that will be merged in
     */
    void mergeMatrices(MatrixZoomDataPP otherMatrixZoom) {
        sum += otherMatrixZoom.sum;
        numRecords += otherMatrixZoom.numRecords;
        for (Map.Entry<Integer, BlockPP> otherBlock : otherMatrixZoom.blocks.entrySet()) {
            int blockNumber = otherBlock.getKey();
            BlockPP block = blocks.get(blockNumber);
            if (block == null) {
                blocks.put(blockNumber, otherBlock.getValue());
                blockNumbers.add(blockNumber);
            } else {
                block.merge(otherBlock.getValue());
                blockNumRecords.put(blockNumber, block.getNumRecords());
            }
        }
        for (int blockNumber : otherMatrixZoom.blockNumbers) {
            blockNumbers.add(blockNumber);
            if (!blockNumRecords.containsKey(blockNumber)) {
                blockNumRecords.put(blockNumber, otherMatrixZoom.blockNumRecords.get(blockNumber));
            } else {
                blockNumRecords.put(blockNumber, blockNumRecords.get(blockNumber) + otherMatrixZoom.blockNumRecords.get(blockNumber));
            }
        }

        tmpFiles.addAll(otherMatrixZoom.tmpFiles);

        for (Map.Entry<Integer, Map<File, Long>> entry : otherMatrixZoom.tmpFilesByBlockNumber.entrySet()) {
            if (!tmpFilesByBlockNumber.containsKey(entry.getKey())) {
                tmpFilesByBlockNumber.put(entry.getKey(), entry.getValue());
            } else {
                for (Map.Entry<File, Long> tmpFile : entry.getValue().entrySet()) {
                    tmpFilesByBlockNumber.get(entry.getKey()).put(tmpFile.getKey(), tmpFile.getValue());
                }
            }
        }
    }

    /**
     * Note -- compressed
     *
     * @param block       Block to write
     * @param sampledData Array to hold a sample of the data (to compute statistics)
     * @throws IOException
     */
    protected void writeBlock(BlockPP block, DownsampledDoubleArrayList sampledData, LittleEndianOutputStream los, Deflater compressor) throws IOException {

        final Map<Point, Float> records = block.getContactRecordMap();//   getContactRecords();

        // System.out.println("Write contact records : records count = " + records.size());

        // Count records first
        int nRecords;
        if (countThreshold > 0) {
            nRecords = 0;
            for (Float value : records.values()) {
                if (value >= countThreshold) {
                    nRecords++;
                }
            }
        } else {
            nRecords = records.size();
        }
        BufferedByteWriter buffer = new BufferedByteWriter(nRecords * 12);
        buffer.putInt(nRecords);
        incrementCellCount(nRecords);

        // Find extents of occupied cells
        int binXOffset = Integer.MAX_VALUE;
        int binYOffset = Integer.MAX_VALUE;
        int binXMax = 0;
        int binYMax = 0;
        for (Map.Entry<Point, Float> entry : records.entrySet()) {
            Point point = entry.getKey();
            binXOffset = Math.min(binXOffset, point.x);
            binYOffset = Math.min(binYOffset, point.y);
            binXMax = Math.max(binXMax, point.x);
            binYMax = Math.max(binYMax, point.y);
        }

        buffer.putInt(binXOffset);
        buffer.putInt(binYOffset);

        // Sort keys in row-major order
        List<Point> keys = new ArrayList<>(records.keySet());
        keys.sort((o1, o2) -> {
            if (o1.y != o2.y) {
                return o1.y - o2.y;
            } else {
                return o1.x - o2.x;
            }
        });
        Point lastPoint = keys.get(keys.size() - 1);
        final short w = (short) (binXMax - binXOffset + 1);
        final int w1 = binXMax - binXOffset + 1;
        final int w2 = binYMax - binYOffset + 1;

        boolean isInteger = true;
        float maxCounts = 0;

        LinkedHashMap<Integer, List<ContactRecord>> rows = new LinkedHashMap<>();
        for (Point point : keys) {
            float counts = records.get(point);
            if (counts >= countThreshold) {
                isInteger = isInteger && (Math.floor(counts) == counts);
                maxCounts = Math.max(counts, maxCounts);

                final int px = point.x - binXOffset;
                final int py = point.y - binYOffset;
                if (!rows.containsKey(py)) {
                    rows.put(py, new ArrayList<>(10));
                }
                List<ContactRecord> row = rows.get(py);
                row.add(new ContactRecord(px, py, counts));
            }
        }

        // Compute size for each representation and choose smallest
        boolean useShort = isInteger && (maxCounts < Short.MAX_VALUE);
        boolean useShortBinX = w1 < Short.MAX_VALUE;
        boolean useShortBinY = w2 < Short.MAX_VALUE;
        int valueSize = useShort ? 2 : 4;

        int lorSize = 0;
        int nDensePts = (lastPoint.y - binYOffset) * w + (lastPoint.x - binXOffset) + 1;


        for (List<ContactRecord> row : rows.values()) {
            lorSize += 4 + row.size() * valueSize;
        }

        buffer.put((byte) (useShort ? 0 : 1));
        buffer.put((byte) (useShortBinX ? 0 : 1));
        buffer.put((byte) (useShortBinY ? 0 : 1));

        //dense calculation is incorrect for v9
        int denseSize = Integer.MAX_VALUE;
        if (lorSize < denseSize) {
            buffer.put((byte) 1);  // List of rows representation
            if (useShortBinY) {
                buffer.putShort((short) rows.size()); // # of rows
            } else {
                buffer.putInt(rows.size());  // # of rows
            }

            for (Map.Entry<Integer, List<ContactRecord>> entry : rows.entrySet()) {

                int py = entry.getKey();
                List<ContactRecord> row = entry.getValue();
                if (useShortBinY) {
                    buffer.putShort((short) py);  // Row number
                } else {
                    buffer.putInt(py); // Row number
                }
                if (useShortBinX) {
                    buffer.putShort((short) row.size());  // size of row
                } else {
                    buffer.putInt(row.size()); // size of row
                }

                for (ContactRecord contactRecord : row) {
                    if (useShortBinX) {
                        buffer.putShort((short) (contactRecord.getBinX()));
                    } else {
                        buffer.putInt(contactRecord.getBinX());
                    }

                    final float counts = contactRecord.getCounts();
                    if (useShort) {
                        buffer.putShort((short) counts);
                    } else {
                        buffer.putFloat(counts);
                    }

                    synchronized(sampledData) {
                        sampledData.add(counts);
                    }
                    incrementSum(counts);
                }
            }

        } else {
            buffer.put((byte) 2);  // Dense matrix

            buffer.putInt(nDensePts);
            buffer.putShort(w);

            int lastIdx = 0;
            for (Point p : keys) {

                int idx = (p.y - binYOffset) * w + (p.x - binXOffset);
                for (int i = lastIdx; i < idx; i++) {
                    // Filler value
                    if (useShort) {
                        buffer.putShort(Short.MIN_VALUE);
                    } else {
                        buffer.putFloat(Float.NaN);
                    }
                }
                float counts = records.get(p);
                if (useShort) {
                    buffer.putShort((short) counts);
                } else {
                    buffer.putFloat(counts);
                }
                lastIdx = idx + 1;

                synchronized(sampledData) {
                    sampledData.add(counts);
                }
                incrementSum(counts);
            }
        }


        byte[] bytes = buffer.getBytes();
        byte[] compressedBytes = compress(bytes, compressor);
        los.write(compressedBytes);

    }

    private synchronized void incrementSum(float counts) {
        sum += counts;
    }

    private synchronized void incrementCellCount(int nRecords) {
        cellCount += nRecords;
    }

    private int addToBlockAndRecordsSets(BlockPP b) {
        int number = b.getNumber();
        blockNumbers.add(number);
        if (blockNumRecords.containsKey(number)) {
            blockNumRecords.put(number, b.getNumRecords() + blockNumRecords.get(number));
        } else {
            blockNumRecords.put(number, b.getNumRecords());
        }
        numRecords += b.getNumRecords();
        return number;
    }
}
