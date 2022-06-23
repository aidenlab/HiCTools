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

package hic.tools.utils.largelists;

import hic.HiCGlobals;
import javastraw.reader.datastructures.ListOfFloatArrays;

import java.util.ArrayList;
import java.util.List;

/**
 * can't use <T> because we need to instantiate the array, otherwise that would have been nice
 */
public class BigFloatsArray {

    final long DEFAULT_LENGTH = BigShortsArray.DEFAULT_LENGTH;
    final long overallLength;
    final List<float[]> internalList = new ArrayList<>();

    public BigFloatsArray(long length) {
        this.overallLength = length;
        long tempLength = length;
        while (tempLength > 0) {
            if (tempLength < DEFAULT_LENGTH) {
                internalList.add(new float[(int) tempLength]);
                break;
            } else {
                internalList.add(new float[(int) DEFAULT_LENGTH]);
                tempLength -= DEFAULT_LENGTH;
            }
        }
    }

    public void clear() {
        internalList.clear();
    }

    public float get(long index) {
        if (index < overallLength) {
            int pseudoRow = (int) (index / DEFAULT_LENGTH);
            int pseudoCol = (int) (index % DEFAULT_LENGTH);
            return internalList.get(pseudoRow)[pseudoCol];
        } else {
            System.err.println("long index exceeds max size of list of arrays while getting: " + index + " " + overallLength);
            Exception ioe = new Exception();
            ioe.printStackTrace();
            return Float.NaN;
        }
    }

    public void set(long index, float value) {
        if (index < overallLength) {
            int pseudoRow = (int) (index / DEFAULT_LENGTH);
            int pseudoCol = (int) (index % DEFAULT_LENGTH);
            internalList.get(pseudoRow)[pseudoCol] = value;
        } else {
            System.err.println("long index exceeds max size of list of arrays while setting");
        }
    }

    public long getLength() {
        return overallLength;
    }

    public BigFloatsArray deepClone() {
        BigFloatsArray clone = new BigFloatsArray(overallLength);
        for (int k = 0; k < internalList.size(); k++) {
            System.arraycopy(internalList.get(k), 0, clone.internalList.get(k), 0, internalList.get(k).length);
        }
        return clone;
    }

    public void multiplyBy(long index, float value) {
        if (index < overallLength) {
            int pseudoRow = (int) (index / DEFAULT_LENGTH);
            int pseudoCol = (int) (index % DEFAULT_LENGTH);
            internalList.get(pseudoRow)[pseudoCol] *= value;
        } else {
            System.err.println("long index exceeds max size of list of arrays while mutiplying");
        }
    }

    public List<float[]> getValues() {
        return internalList;
    }

    public ListOfFloatArrays convertToRegular() {
        ListOfFloatArrays clone = new ListOfFloatArrays(overallLength);
        for (long k = 0; k < getLength(); k++) {
            clone.set(k, get(k));
        }
        return clone;
    }

    private int getNumThreads() {
        return Math.min(HiCGlobals.numCPUMatrixThreads, internalList.size());
    }

    public BigDoublesArray convertToDoubles() {
        BigDoublesArray clone = new BigDoublesArray(overallLength);
        for (int k = 0; k < internalList.size(); k++) {
            double[] dest = clone.internalList.get(k);
            float[] src = internalList.get(k);
            for (int q = 0; q < dest.length; q++) {
                dest[q] = src[q];
            }
        }
        return clone;
    }
}

