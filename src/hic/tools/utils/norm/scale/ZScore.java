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

package hic.tools.utils.norm.scale;

import javastraw.reader.datastructures.ListOfIntArrays;

public class ZScore {
    private final static int UPPER_CUTOFF = 2;
    private final float mu, std;

    public ZScore(ListOfIntArrays array) {
        float initialMu = getMean(array);
        float initialStd = getStandardDev(array, initialMu);
        mu = getConservativeMean(array, initialMu, initialStd);
        std = getConservativeStd(array, mu, initialMu, initialStd);
    }

    private float getConservativeStd(ListOfIntArrays array, float mu, float initialMu, float initialStd) {
        double std = 0;
        long num = 0;
        for (long p = 0; p < array.getLength(); p++) {
            int valP = array.get(p);
            if (valP > 0 && getZ(valP, initialMu, initialStd) < UPPER_CUTOFF) {
                float diff = valP - mu;
                std += diff * diff;
                num++;
            }
        }

        return (float) Math.sqrt(std / Math.max(num, 1));
    }

    private float getConservativeMean(ListOfIntArrays array, float initialMu, float initialStd) {
        double mu = 0;
        long num = 0;
        for (long p = 0; p < array.getLength(); p++) {
            int valP = array.get(p);
            if (valP > 0 && getZ(valP, initialMu, initialStd) < UPPER_CUTOFF) {
                mu += valP;
                num++;
            }
        }
        return (float) (mu / num);
    }

    private static float getMean(ListOfIntArrays array) {
        double mu = 0;
        long num = 0;
        for (long p = 0; p < array.getLength(); p++) {
            int valP = array.get(p);
            if (valP > 0) {
                mu += valP;
                num++;
            }
        }
        return (float) (mu / num);
    }

    private static float getStandardDev(ListOfIntArrays array, float mu) {
        double std = 0;
        long num = 0;
        for (long p = 0; p < array.getLength(); p++) {
            int valP = array.get(p);
            if (valP > 0) {
                float diff = valP - mu;
                std += diff * diff;
                num++;
            }
        }

        return (float) Math.sqrt(std / Math.max(num, 1));
    }

    public float getCutoff(float z) {
        return mu + z * std;
    }

    private float getZ(float val, float mu, float std) {
        return (val - mu) / std;
    }
}
