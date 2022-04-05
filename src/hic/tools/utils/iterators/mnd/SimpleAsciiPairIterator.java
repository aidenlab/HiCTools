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

package hic.tools.utils.iterators.mnd;

import hic.HiCGlobals;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.zip.GZIPInputStream;

public class SimpleAsciiPairIterator extends GenericPairIterator implements PairIterator {

    public SimpleAsciiPairIterator(String path) throws IOException {
        super(new MNDFileParser(new SimpleLineParser()));
        if (path.endsWith(".gz")) {
            InputStream gzipStream = new GZIPInputStream(new FileInputStream(path));
            Reader decoder = new InputStreamReader(gzipStream, StandardCharsets.UTF_8);
            this.reader = new BufferedReader(decoder, 4194304);
        } else {
            this.reader = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(path)), HiCGlobals.bufferSize);
        }
        advance();
    }

    public String getChromosomeNameFromIndex(int chrIndex) {
        return mndFileParser.getChromosomeNameFromIndex(chrIndex);
    }
}
