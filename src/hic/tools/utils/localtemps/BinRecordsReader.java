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

import htsjdk.tribble.util.LittleEndianInputStream;
import javastraw.reader.block.ContactRecord;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Iterator;

public class BinRecordsReader implements Iterator<ContactRecord> {

    protected final LittleEndianInputStream is;
    protected ContactRecord next;

    public BinRecordsReader(String path) throws IOException {
        is = new LittleEndianInputStream(new BufferedInputStream(Files.newInputStream(Paths.get(path))));
        advance();
    }

    public boolean hasNext() {
        return next != null;
    }

    public ContactRecord next() {
        ContactRecord retValue = next;
        advance();
        return retValue;
    }

    public void remove() {
    }

    public void close() {
        if (is != null) try {
            is.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    protected void advance() {
        try {
            int binX = is.readInt();
            int binY = is.readInt();
            float val = is.readFloat();
            next = new ContactRecord(binX, binY, val);
        } catch (IOException e) {
            next = null;
            if (!(e instanceof EOFException)) {
                e.printStackTrace();
            }
        }
    }
}

