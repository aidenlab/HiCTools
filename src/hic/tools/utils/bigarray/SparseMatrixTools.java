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
import hic.tools.utils.largelists.BigIntsArray;
import hic.tools.utils.original.ExpectedValueCalculation;
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

    public static void matrixVectorMult(BigIntsArray vector,
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

    public static void populateNormedExpected(int chrIdx, ListOfFloatArrays expectedVector, ExpectedValueCalculation exp, int x, int y, float counts) {
        if (counts > 0) {
            final float vx = expectedVector.get(x);
            final float vy = expectedVector.get(y);
            if (vx > 0 && vy > 0) {
                double value = counts / (vx * vy);
                if (value > 0) {
                    exp.addDistance(chrIdx, x, y, value);
                }
            }
        }
    }

    public static void sumRawAndNorm(double[] normSum, double[] sum, int x, int y, float counts,
                                     ListOfFloatArrays newNormVector) {
        double valX = newNormVector.get(x);
        double valY = newNormVector.get(y);

        if (valX > 0 && valY > 0) {
            double normalizedValue = counts / (valX * valY);
            normSum[0] += normalizedValue;
            sum[0] += counts;
            if (x != y) {
                normSum[0] += normalizedValue;
                sum[0] += counts;
            }
        }
    }

    public static void sumScaleFactor(ListOfFloatArrays norm, double[] mSum, double[] nSum, int x, int y, float value) {
        double valX = norm.get(x);
        double valY = norm.get(y);
        if (valX > 0 && valY > 0) {
            // want total sum of matrix, not just upper triangle
            if (x == y) {
                nSum[0] += value / (valX * valY);
                mSum[0] += value;
            } else {
                nSum[0] += 2 * value / (valX * valY);
                mSum[0] += 2 * value;
            }
        }
    }

    public static void updateRowSums(ListOfFloatArrays sums, int x, int y, float value) {
        sums.addTo(x, value);
        if (x != y) {
            sums.addTo(y, value);
        }
    }
}
