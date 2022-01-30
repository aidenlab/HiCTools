/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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

package juicebox.tools.utils.norm;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.iterators.IteratorContainer;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import juicebox.tools.utils.norm.scale.ScaleHandler;

public class CustomNormalizationVector extends NormalizationVector {
    private final boolean needsToBeScaledTo;

    public CustomNormalizationVector(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit, int resolution,
                                     ListOfDoubleArrays data, boolean needsToBeScaledTo) {
        super(type, chrIdx, unit, resolution, data);
        this.needsToBeScaledTo = needsToBeScaledTo;
    }

    public boolean doesItNeedToBeScaledTo() {
        return needsToBeScaledTo;
    }

    public NormalizationVector mmbaScaleToVector(Dataset ds) {
        Chromosome chromosome = ds.getChromosomeHandler().getChromosomeFromIndex(chrIdx);
        MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chromosome, chromosome, new HiCZoom(unit, resolution));
        if (zd == null) return null;
        return mmbaScaleToVector(zd.getIteratorContainer());
    }

    public NormalizationVector mmbaScaleToVector(IteratorContainer ic) {
        ListOfFloatArrays newNormVector = ScaleHandler.scale(ic, data.convertToFloats(), resolution);
        ScaleHandler.normalizeVectorByScaleFactor(newNormVector, ic);
        ListOfDoubleArrays newDoubleNormVector = newNormVector.convertToDoubles();
        return new NormalizationVector(type, chrIdx, unit, resolution, newDoubleNormVector);
    }
}
