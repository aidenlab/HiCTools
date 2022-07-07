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

package hic.tools.utils.norm;

import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

public class FloatNormVector {
    private final NormalizationType type;
    private final int chrIdx;
    private final HiCZoom.HiCUnit unit;
    private final int resolution;
    private final ListOfFloatArrays data;

    public FloatNormVector(NormalizationType type, int chrIdx, HiCZoom zoom, ListOfFloatArrays data) {
        this.type = type;
        this.chrIdx = chrIdx;
        this.unit = zoom.getUnit();
        this.resolution = zoom.getBinSize();
        this.data = data;
    }

    public static String getKey(NormalizationType type, int chrIdx, String unit, int resolution) {
        return type + "_" + chrIdx + "_" + unit + "_" + resolution;
    }

    public int getChrIdx() {
        return this.chrIdx;
    }

    public int getResolution() {
        return this.resolution;
    }

    public String getKey() {
        return getKey(this.type, this.chrIdx, this.unit.toString(), this.resolution);
    }

    public ListOfFloatArrays getData() {
        return this.data;
    }

    public NormalizationType getNormType() {
        return type;
    }
}