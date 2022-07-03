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

package hic.tools.utils.localtemps;

import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.iterators.GenomeWideIterator;
import javastraw.reader.type.HiCZoom;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class BinRecordsWriter {

    private static int internalCount = 0;

    public static void saveAllContacts(Iterator<ContactRecord> iterator, int limit,
                                       List<String> filenames, String stem) throws IOException {
        LittleEndianOutputStream les = createNewTempFile(filenames, stem);
        int counter = 0;
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            les.writeInt(record.getBinX());
            les.writeInt(record.getBinY());
            les.writeFloat(record.getCounts());
            counter++;
            if (counter >= limit) {
                les.close();
                les = createNewTempFile(filenames, stem);
                counter = 0;
            }
        }
        if (counter > 0) {
            les.close();
        }
    }

    private static LittleEndianOutputStream createNewTempFile(List<String> files, String stem) throws IOException {
        File tempFile = File.createTempFile(stem + (internalCount++), ".tmp.bin");
        tempFile.deleteOnExit();
        files.add(tempFile.getAbsolutePath());

        BufferedOutputStream bos = new BufferedOutputStream(Files.newOutputStream(tempFile.toPath()));
        return new LittleEndianOutputStream(bos);
    }

    public static void saveAllGWContacts(Dataset ds, ChromosomeHandler handler, HiCZoom zoom, boolean includeIntra,
                                         int limit, List<String> filenames,
                                         String intraStem, String interStem) {

        ExecutorService executor = Executors.newFixedThreadPool(2);

        Runnable worker = () -> {
            try {
                saveAllContacts(new GenomeWideIterator(ds, handler, zoom, false, true),
                        limit, filenames, interStem);
            } catch (IOException e) {
                System.err.println("Unable to save inter data to local files");
                e.printStackTrace();
                System.exit(31);
            }
        };
        executor.execute(worker);

        if (includeIntra) {
            Runnable worker2 = () -> {
                try {
                    saveAllContacts(new GenomeWideIterator(ds, handler, zoom, true, false),
                            limit, filenames, intraStem);
                } catch (IOException e) {
                    System.err.println("Unable to save intra data to local files");
                    e.printStackTrace();
                    System.exit(31);
                }
            };
            executor.execute(worker2);
        }

        executor.shutdown();

        while (!executor.isTerminated()) {
        }
    }
}
