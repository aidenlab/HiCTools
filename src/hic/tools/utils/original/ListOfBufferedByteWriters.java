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
import org.broad.igv.tdf.BufferedByteWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ListOfBufferedByteWriters {

    List<BufferedByteWriter> bufferList = new ArrayList<>();

    ListOfBufferedByteWriters() {
        expandBuffer();
    }

    public void putInt(int value) throws IOException {
        bufferList.get(index()).putInt(value);
    }

    public void putLong(long value) throws IOException {
        bufferList.get(index()).putLong(value);
    }

    public void putFloat(float value) throws IOException {
        bufferList.get(index()).putFloat(value);
    }

    public void putNullTerminatedString(String value) throws IOException {
        bufferList.get(index()).putNullTerminatedString(value);
    }

    public void expandBufferIfNeeded(int buffer) {
        if (Integer.MAX_VALUE - bufferList.get(index()).bytesWritten() < buffer) {
            bufferList.add(new BufferedByteWriter());
        }
    }

    private int index() {
        return bufferList.size() - 1;
    }

    public void expandBuffer() {
        bufferList.add(new BufferedByteWriter());
    }

    public long getTotalBytes() {
        long total = 0;
        for (BufferedByteWriter writer : bufferList) {
            total += writer.getBytes().length;
        }
        return total;
    }

    public void writeToOutput(LittleEndianOutputStream los) throws IOException {
        for (BufferedByteWriter writer : bufferList) {
            los.write(writer.getBytes());
        }
    }
}
