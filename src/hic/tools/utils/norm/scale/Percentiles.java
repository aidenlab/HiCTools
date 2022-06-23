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

package hic.tools.utils.norm.scale;

import javastraw.reader.datastructures.ListOfIntArrays;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.List;

public class Percentiles {
    private final double[] cutoffs;

    public Percentiles(ListOfIntArrays array, float startingPercentile, float deltaPercentile, float maxPercentile) {
        cutoffs = createCutoffs(getStats(array), startingPercentile, deltaPercentile, maxPercentile);
    }

    private static DescriptiveStatistics getStats(ListOfIntArrays array) {
        DescriptiveStatistics statistics = new DescriptiveStatistics();
        for (long p = 0; p < array.getLength(); p++) {
            int valP = array.get(p);
            if (valP > 0) {
                statistics.addValue(valP);
            }
        }
        return statistics;
    }

    private double[] createCutoffs(DescriptiveStatistics stats, float startingPercentile, float deltaPercentile, float maxPercentile) {
        List<Double> values = new ArrayList<>();
        for (float cutoff = startingPercentile; cutoff < maxPercentile; cutoff += deltaPercentile) {
            values.add(stats.getPercentile(cutoff));
        }
        return toArray(values);
    }

    private double[] toArray(List<Double> values) {
        double[] array = new double[values.size()];
        for (int i = 0; i < array.length; i++) {
            array[i] = values.get(i);
        }
        return array;
    }

    public double get(int cutoffIndex) {
        if (cutoffIndex < cutoffs.length) {
            return cutoffs[cutoffIndex];
        }
        return cutoffs[cutoffs.length - 1];
    }

    public int getMaxNumCutoffs() {
        return cutoffs.length;
    }
}
