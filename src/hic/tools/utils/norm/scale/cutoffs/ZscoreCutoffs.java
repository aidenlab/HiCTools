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

package hic.tools.utils.norm.scale.cutoffs;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.reader.datastructures.ListOfIntArrays;

import java.util.ArrayList;
import java.util.List;

public class ZscoreCutoffs extends Cutoffs {

    public ZscoreCutoffs(ListOfIntArrays array, float startingZscore, float deltaZscore, float maxZscore) {
        super(createCutoffs(getStats(array), startingZscore, deltaZscore, maxZscore));
    }

    private static Zscore getStats(ListOfIntArrays array) {
        Welford welford = new Welford();
        for (long p = 0; p < array.getLength(); p++) {
            int valP = array.get(p);
            if (valP > 0) {
                welford.addValue(valP);
            }
        }
        return welford.getZscore();
    }

    private static double[] createCutoffs(Zscore zscore, float minZ, float deltaZ, float maxZ) {
        List<Double> values = new ArrayList<>();
        values.add(0.0);
        for (float cutoff = minZ; cutoff < maxZ; cutoff += deltaZ) {
            values.add(zscore.getValForZscore(cutoff));
        }
        return toArray(values);
    }
}
