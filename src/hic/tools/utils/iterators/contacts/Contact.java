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

package hic.tools.utils.iterators.contacts;

import javastraw.reader.block.ContactRecord;

public class Contact {

    private final int chr1;
    private final int pos1;
    private final int chr2;
    private final int pos2;
    private final float score;  // The score (or count)

    public Contact(int chr1Index, int chr2Index, ContactRecord record, int resolution) {
        this.chr1 = chr1Index;
        this.chr2 = chr2Index;
        this.pos1 = record.getBinX() * resolution;
        this.pos2 = record.getBinY() * resolution;
        this.score = record.getCounts();
    }


    public int getChr1() {
        return chr1;
    }

    public int getPos1() {
        return pos1;
    }

    public int getChr2() {
        return chr2;
    }

    public int getPos2() {
        return pos2;
    }

    public float getScore() {
        return score;
    }
}
