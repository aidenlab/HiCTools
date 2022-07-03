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

package hic.tools.utils.bigarray;

import hic.tools.utils.largelists.BigDoublesArray;
import hic.tools.utils.largelists.BigFloatsArray;
import hic.tools.utils.largelists.BigShortsArray;
import javastraw.reader.datastructures.ListOfFloatArrays;

public class SparseMatrixTools {
    public static void matrixVectorMult(BigFloatsArray vector,
                                        BigDoublesArray sumVector, int x, int y, float c) {
        double counts = c;
        if (x == y) {
            counts *= .5;
        }
        sumVector.addTo(x, counts * vector.get(y));
        sumVector.addTo(y, counts * vector.get(x));
    }

    public static void matrixVectorMult(BigShortsArray vector,
                                        BigDoublesArray sumVector, int x, int y, float c) {
        double counts = c;
        if (x == y) {
            counts *= .5;
        }
        sumVector.addTo(x, counts * vector.get(y));
        sumVector.addTo(y, counts * vector.get(x));
    }

    public static void invertVector(ListOfFloatArrays newNormVector) {
        for (long k = 0; k < newNormVector.getLength(); k++) {
            float kVal = newNormVector.get(k);
            if (kVal > 0) {
                newNormVector.set(k, 1.f / kVal);
            } else {
                newNormVector.set(k, Float.NaN);
            }
        }
    }
}
