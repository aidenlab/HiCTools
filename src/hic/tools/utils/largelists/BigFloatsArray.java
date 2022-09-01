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

import com.google.common.util.concurrent.AtomicDouble;
import hic.HiCGlobals;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.tools.ParallelizationTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

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

    public void add(long index, float value) {
        if (index < overallLength) {
            int pseudoRow = (int) (index / DEFAULT_LENGTH);
            int pseudoCol = (int) (index % DEFAULT_LENGTH);
            internalList.get(pseudoRow)[pseudoCol] += value;
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

    public static double parCalculateError(BigFloatsArray col, BigFloatsArray scale, BigShortsArray target, BigShortsArray bad) {
        AtomicDouble atomicDouble = new AtomicDouble(0);
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(col.getNumThreads(), () -> {
            int i = index.getAndIncrement();
            float err = 0;
            while (i < col.internalList.size()) {
                float[] colA = col.internalList.get(i);
                float[] scaleA = scale.internalList.get(i);
                short[] targetA = target.internalList.get(i);
                short[] badA = bad.internalList.get(i);

                for (int z = 0; z < colA.length; z++) {
                    if (badA[z] == 1) continue;
                    float tempErr = Math.abs((colA[z] * scaleA[z] - targetA[z]));
                    if (tempErr > err) {
                        err = tempErr;
                    }
                }
                i = index.getAndIncrement();
            }
            synchronized (atomicDouble) {
                if (err > atomicDouble.get()) {
                    atomicDouble.set(err);
                }
            }
        });
        return atomicDouble.get();
    }

    public static double calculateError90(BigFloatsArray col, BigFloatsArray scale,
                                          BigShortsArray target, BigShortsArray bad) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < col.internalList.size(); i++) {
            float[] colA = col.internalList.get(i);
            float[] scaleA = scale.internalList.get(i);
            short[] targetA = target.internalList.get(i);
            short[] badA = bad.internalList.get(i);

            for (int z = 0; z < colA.length; z++) {
                if (badA[z] == 1) continue;
                float tempErr = Math.abs((colA[z] * scaleA[z] - targetA[z]));
                stats.addValue(tempErr);
            }
        }
        return stats.getPercentile(90);
    }


    public static double parCalculateConvergenceError(BigFloatsArray calculatedVectorB, BigFloatsArray current,
                                                      BigShortsArray bad) {
        AtomicDouble atomicDouble = new AtomicDouble(0);
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(current.getNumThreads(), () -> {
            int i = index.getAndIncrement();
            double err = 0;
            while (i < current.internalList.size()) {

                float[] calcA = calculatedVectorB.internalList.get(i);
                float[] currA = current.internalList.get(i);
                short[] badA = bad.internalList.get(i);

                for (int z = 0; z < calcA.length; z++) {
                    if (badA[z] == 1) continue;
                    double relativeErr = Math.abs(((calcA[z] - currA[z]) / (calcA[z] + currA[z])));
                    //float tempErr = Math.abs(calcA[z] - currA[z]);
                    if (relativeErr > err) {
                        err = relativeErr;
                    }
                }
                i = index.getAndIncrement();
            }
            synchronized (atomicDouble) {
                if (err > atomicDouble.get()) {
                    atomicDouble.set(err);
                }
            }
        });
        return atomicDouble.get();
    }

    public void parSetToGeoMean(BigFloatsArray a, BigFloatsArray b) {
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int i = index.getAndIncrement();
            while (i < internalList.size()) {
                float[] result = internalList.get(i);
                float[] a1 = a.internalList.get(i);
                float[] b1 = b.internalList.get(i);
                for (int p = 0; p < result.length; p++) {
                    result[p] = (float) Math.sqrt(a1[p] * b1[p]);
                }
                i = index.getAndIncrement();
            }
        });
    }

    public void parSetTo(BigFloatsArray srcArrays) {
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int i = index.getAndIncrement();
            while (i < internalList.size()) {
                float[] dest = internalList.get(i);
                float[] src = srcArrays.internalList.get(i);
                System.arraycopy(src, 0, dest, 0, dest.length);
                i = index.getAndIncrement();
            }
        });
    }

    public void parMultiplyByOneMinus(BigShortsArray array) {
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int i = index.getAndIncrement();
            while (i < internalList.size()) {
                float[] orig = internalList.get(i);
                short[] arr = array.internalList.get(i);
                for (int p = 0; p < orig.length; p++) {
                    orig[p] *= (1 - arr[p]);
                }
                i = index.getAndIncrement();
            }
        });
    }

    public void parMultiplyBy(BigFloatsArray dv) {
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int i = index.getAndIncrement();
            while (i < internalList.size()) {
                float[] orig = internalList.get(i);
                float[] arr = dv.internalList.get(i);
                for (int p = 0; p < orig.length; p++) {
                    orig[p] *= arr[p];
                }
                i = index.getAndIncrement();
            }
        });
    }

    public void parSetToDivision(BigShortsArray num, BigFloatsArray denom) {
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int i = index.getAndIncrement();
            while (i < internalList.size()) {
                float[] orig = internalList.get(i);
                short[] num1 = num.internalList.get(i);
                float[] denom1 = denom.internalList.get(i);
                for (int p = 0; p < orig.length; p++) {
                    orig[p] = num1[p] / denom1[p];
                }
                i = index.getAndIncrement();
            }
        });
    }

    private int getNumThreads() {
        return Math.min(HiCGlobals.normThreads, internalList.size());
    }

    public void parScaleByRatio(BigShortsArray num, BigFloatsArray denom) {
        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
            int i = index.getAndIncrement();
            while (i < internalList.size()) {
                float[] orig = internalList.get(i);
                short[] num1 = num.internalList.get(i);
                float[] denom1 = denom.internalList.get(i);
                for (int p = 0; p < orig.length; p++) {
                    orig[p] *= (num1[p] / denom1[p]);
                }
                i = index.getAndIncrement();
            }
        });
    }

    public void addValuesFrom(BigFloatsArray other) {
        if (overallLength == other.overallLength) {
            for (int i = 0; i < internalList.size(); i++) {
                for (int j = 0; j < internalList.get(i).length; j++) {
                    internalList.get(i)[j] += other.internalList.get(i)[j];
                }
            }
        } else {
            System.err.println("Adding objects of different sizes!");
        }
    }
}

