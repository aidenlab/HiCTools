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

package hic.tools.utils.norm;

import hic.tools.utils.largelists.BigListOfByteWriters;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.expected.ExpectedValueFunctionImpl;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import org.broad.igv.tdf.BufferedByteWriter;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class NormVectorUpdater {

    static void updateNormVectorIndexWithVector(List<NormalizationVectorIndexEntry> normVectorIndex,
                                                BigListOfByteWriters bufferList, ListOfFloatArrays vec,
                                                int chrIdx, NormalizationType type, HiCZoom zoom) throws IOException {
        long position = bufferList.getBytesWritten();

        putFloatArraysIntoBufferList(bufferList, vec.getValues());

        long newPos = bufferList.getBytesWritten();
        int sizeInBytes = (int) (newPos - position);
        normVectorIndex.add(new NormalizationVectorIndexEntry(type.toString(), chrIdx, zoom.getUnit().toString(),
                zoom.getBinSize(), position, sizeInBytes));
    }

    static void putFloatArraysIntoBufferList(BigListOfByteWriters bufferList,
                                             List<float[]> arrays) throws IOException {

        bufferList.expandBufferIfNeeded(8);
        long vectorLength = 0;
        for (float[] array : arrays) {
            vectorLength += array.length;
        }
        bufferList.putLong(vectorLength);

        for (float[] array : arrays) {
            bufferList.expandBufferIfNeeded(4 * array.length);
            for (float val : array) {
                bufferList.putFloat(val);
            }
        }
    }

    private static void putMapValuesIntoBuffer(BigListOfByteWriters bufferList, Map<Integer, Double> hashmap) throws IOException {
        int bytesNeeded = 4 + (8 * hashmap.size());
        bufferList.expandBufferIfNeeded(bytesNeeded);
        bufferList.putInt(hashmap.size());
        List<Integer> keys = new ArrayList<>(hashmap.keySet());
        Collections.sort(keys);
        for (Integer key : keys) {
            bufferList.putInt(key);
            bufferList.putFloat(hashmap.get(key).floatValue());
        }
    }

    private static void writeExpectedToBuffer(RandomAccessFile raf,
                                              BigListOfByteWriters bufferList, long filePosition) throws IOException {
        raf.getChannel().position(filePosition);
        bufferList.writeToRAF(raf);
    }

    private static void handleVersionSix(RandomAccessFile raf, int version) throws IOException {
        if (version < 6) {
            // Update version
            // Master index
            raf.getChannel().position(4);
            BufferedByteWriter buffer = new BufferedByteWriter();
            buffer.putInt(6);
            raf.write(buffer.getBytes());
        }
    }

    /**
     * Compute the size of the index in bytes.  This is needed to set offsets for the actual index entries.  The
     * easiest way to do this is to write it to a buffer and check the size
     *
     * @param buffer          Buffer to write to
     * @param normVectorIndex Normalization index to write
     */
    static void writeNormIndex(BufferedByteWriter buffer, List<NormalizationVectorIndexEntry> normVectorIndex) throws IOException {
        buffer.putInt(normVectorIndex.size());
        for (NormalizationVectorIndexEntry entry : normVectorIndex) {
            buffer.putNullTerminatedString(entry.type);
            buffer.putInt(entry.chrIdx);
            buffer.putNullTerminatedString(entry.unit);
            buffer.putInt(entry.resolution);
            buffer.putLong(entry.position);
            buffer.putLong(entry.sizeInBytes);
        }
    }

    static void writeNormsToUpdateFile(DatasetReaderV2 reader, String path, boolean useCalcNotFunc,
                                       List<ExpectedValueCalculation> expectedValueCalculations,
                                       Map<String, ExpectedValueFunction> expectedValueFunctionMap,
                                       List<NormalizationVectorIndexEntry> normVectorIndices,
                                       BigListOfByteWriters normVectorBuffers, String message) throws IOException {
        int version = reader.getVersion();
        long filePosition = reader.getNormFilePosition();
        long nviHeaderPosition = reader.getNviHeaderPosition();


        try (RandomAccessFile raf = new RandomAccessFile(path, "rw")) {
            handleVersionSix(raf, version);
            BigListOfByteWriters bufferList = new BigListOfByteWriters();

            if (useCalcNotFunc) {
                writeExpectedValues(bufferList, expectedValueCalculations);
            } else {
                writeExpectedValues(bufferList, expectedValueFunctionMap);
            }

            writeExpectedToBuffer(raf, bufferList, filePosition);
            writeNormsToBuffer(raf, normVectorIndices, normVectorBuffers, nviHeaderPosition);
        }

        System.out.println(message);
    }

    private static void writeExpectedValues(BigListOfByteWriters bufferList,
                                            List<ExpectedValueCalculation> expectedValueCalculations) throws IOException {
        bufferList.expandBufferIfNeeded(4);
        bufferList.putInt(expectedValueCalculations.size());

        for (ExpectedValueCalculation ev : expectedValueCalculations) {
            ev.computeDensity();
            HiCZoom.HiCUnit unit = HiCZoom.HiCUnit.BP;
            appendExpectedValuesToBuffer(bufferList, ev.getType(),
                    unit, ev.getGridSize(), ev.getDensityAvg(),
                    ev.getChrScaleFactors());
        }
    }

    private static void writeExpectedValues(BigListOfByteWriters bufferList, Map<String, ExpectedValueFunction> expectedValueFunctionMap) throws IOException {

        bufferList.expandBufferIfNeeded(4);
        bufferList.putInt(expectedValueFunctionMap.size());

        for (ExpectedValueFunction function : expectedValueFunctionMap.values()) {
            appendExpectedValuesToBuffer(bufferList, function.getNormalizationType(),
                    function.getUnit(), function.getBinSize(),
                    function.getExpectedValuesNoNormalization(),
                    ((ExpectedValueFunctionImpl) function).getNormFactors());
        }
    }

    private static void appendExpectedValuesToBuffer(BigListOfByteWriters bufferList,
                                                     NormalizationType normalizationType,
                                                     HiCZoom.HiCUnit unit, int binSize,
                                                     ListOfDoubleArrays expectedValuesNoNormalization,
                                                     Map<Integer, Double> normFactors) throws IOException {

        int bytesNeeded = normalizationType.toString().length() + 1;
        bytesNeeded += unit.toString().length() + 1;
        bytesNeeded += 4;

        bufferList.expandBufferIfNeeded(bytesNeeded);
        bufferList.putNullTerminatedString(normalizationType.toString());
        bufferList.putNullTerminatedString(unit.toString());
        bufferList.putInt(binSize);

        putFloatArraysIntoBufferList(bufferList,
                expectedValuesNoNormalization.convertToFloats().getValues());

        putMapValuesIntoBuffer(bufferList, normFactors);
    }

    private static void writeNormsToBuffer(RandomAccessFile raf, List<NormalizationVectorIndexEntry> normVectorIndex,
                                           BigListOfByteWriters normVectorBuffers, long nviHeaderPosition) throws IOException {
        // Get the size of the index in bytes, to compute an offset for the actual entries.
        BufferedByteWriter buffer = new BufferedByteWriter();
        writeNormIndex(buffer, normVectorIndex);
        long normVectorStartPosition = raf.getChannel().position() + buffer.bytesWritten();
        long size = buffer.bytesWritten();
        long NVI = normVectorStartPosition - size;
        // write NVI, size
        raf.getChannel().position(nviHeaderPosition);

        BufferedByteWriter headerBuffer = new BufferedByteWriter();
        headerBuffer.putLong(NVI);
        headerBuffer.putLong(size);
        raf.write(headerBuffer.getBytes());

        // reset pointer to where we were
        raf.getChannel().position(NVI);

        // Update index entries
        for (NormalizationVectorIndexEntry entry : normVectorIndex) {
            entry.position += normVectorStartPosition;
        }

        // Now write for real
        buffer = new BufferedByteWriter();
        writeNormIndex(buffer, normVectorIndex);
        raf.write(buffer.getBytes());
        // Finally the norm vectors
        normVectorBuffers.writeToRAF(raf);
    }
}
