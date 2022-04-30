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

import java.awt.*;
import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.Deflater;

public class RecordBlockUtils {

    public static int getNumberOfRecords(Map<Point, Float> records, int countThreshold) {
        if (countThreshold > 0) {
            int nRecords = 0;
            for (Float value : records.values()) {
                if (value >= countThreshold) {
                    nRecords++;
                }
            }
            return nRecords;
        }
        return records.size();
    }

    public static byte[] compress(byte[] data, Deflater compressor) {
        compressor.reset();
        compressor.setInput(data);
        compressor.finish();

        // Create an expandable byte array to hold the compressed data.
        // You cannot use an array that's the same size as the original because
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

    public static void readAndMerge(BlockPP currentBlock, Map.Entry<File, Long> entry) throws IOException {
        BlockPP tmpBlock = readTmpBlock(entry.getKey(), entry.getValue());
        if (tmpBlock != null) {
            currentBlock.merge(tmpBlock);
        }
    }

    public static BlockPP readTmpBlock(File file, long filePosition) throws IOException {
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

    public static Map<Point, Float> readContactRecordsToMap(int nRecords, byte[] bytes) throws IOException {
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
}
