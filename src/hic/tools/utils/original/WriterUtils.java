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

import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.IOException;
import java.util.zip.Deflater;

public class WriterUtils {
    public static Deflater getDefaultCompressor() {
        Deflater compressor = new Deflater();
        compressor.setLevel(Deflater.DEFAULT_COMPRESSION);
        return compressor;
    }

    public static void writeZoomHeader(MatrixZoomDataPP zd, LittleEndianOutputStream los) throws IOException {
        int numberOfBlocks = zd.blockNumbers.size();
        los.writeString(zd.getUnit().toString());  // Unit
        los.writeInt(zd.getZoom());     // zoom index,  lowest res is zero
        los.writeFloat((float) zd.getSum());      // sum
        los.writeFloat((float) zd.getOccupiedCellCount());
        los.writeFloat((float) zd.getPercent5());
        los.writeFloat((float) zd.getPercent95());
        los.writeInt(zd.getBinSize());
        los.writeInt(zd.getBlockBinCount());
        los.writeInt(zd.getBlockColumnCount());
        los.writeInt(numberOfBlocks);
        zd.blockIndexPosition = los.getWrittenCount();
        // Placeholder for block index
        for (int i = 0; i < numberOfBlocks; i++) {
            los.writeInt(0);
            los.writeLong(0L);
            los.writeInt(0);
        }
    }
}
