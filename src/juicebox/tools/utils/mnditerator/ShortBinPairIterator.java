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

package juicebox.tools.utils.mnditerator;

import java.io.EOFException;
import java.io.IOException;

public class ShortBinPairIterator extends BinPairIterator {

    public ShortBinPairIterator(String path) throws IOException {
        super(path);
    }

    @Override
    protected void advance() {
        try {
            int chr1 = is.readInt();
            int pos1 = is.readInt();
            int chr2 = is.readInt();
            int pos2 = is.readInt();
            next = new AlignmentPair(chr1, pos1, chr2, pos2);

            float score = is.readFloat();
            next.setScore(score);
        } catch (IOException e) {
            next = null;
            if (!(e instanceof EOFException)) {
                e.printStackTrace();
            }
        }
    }
}
