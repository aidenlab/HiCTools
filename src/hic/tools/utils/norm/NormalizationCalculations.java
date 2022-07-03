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

import hic.tools.clt.old.NormalizationBuilder;
import hic.tools.utils.bigarray.BigContactList;
import hic.tools.utils.largelists.BigFloatsArray;
import hic.tools.utils.norm.scale.FinalScale;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.type.NormalizationType;

public class NormalizationCalculations {

    private final long matrixSize; // x and y symmetric
    //private boolean isEnoughMemory = false;
    private final BigContactList ba;
    private final int resolution;

    public NormalizationCalculations(BigContactList ba, int resolution) {
        this.ba = ba;
        this.matrixSize = ba.getMatrixSize();
        this.resolution = resolution;
    }

    public ListOfFloatArrays getNorm(NormalizationType normOption, String stem) {
        ListOfFloatArrays norm;
        if (NormalizationBuilder.usesVC(normOption)) {
            norm = computeVC();
        } else if (NormalizationBuilder.usesSCALE(normOption)) {
            norm = computeSCALE(computeVC(), stem);
        } else if (NormalizationBuilder.isNONE(normOption)) {
            return new ListOfFloatArrays(matrixSize, 1);
        } else {
            System.err.println("Not supported for normalization " + normOption);
            return null;
        }

        if (norm != null && norm.getLength() > 0) {
            double factor = getSumFactor(norm);
            norm.multiplyEverythingBy(factor);
        }
        return norm;
    }
    
    /**
     * Compute vanilla coverage norm, just the sum of the rows
     *
     * @return Normalization vector
     */
    ListOfFloatArrays computeVC() {
        return ba.getRowSums();
    }

    /**
     * Get the sum of the normalized matrix
     *
     * @param norm Normalization vector
     * @return Square root of ratio of original to normalized vector
     */
    public double getSumFactor(ListOfFloatArrays norm) {
        double[] normMatrixSums = getNormMatrixSumFactor(norm);
        return Math.sqrt(normMatrixSums[0] / normMatrixSums[1]);
    }
    
    public double[] getNormMatrixSumFactor(ListOfFloatArrays norm) {
        return ba.getNormMatrixSumFactor(norm);
    }

    private BigFloatsArray getInitialStartingVector(ListOfFloatArrays vc) {
        BigFloatsArray initial = new BigFloatsArray(vc.getLength());
        for (long i = 0; i < vc.getLength(); i++) {
            initial.set(i, (float) Math.sqrt(vc.get(i)));
        }
        return initial;
    }

    public ListOfFloatArrays computeSCALE(ListOfFloatArrays vc, String stem) {
        BigFloatsArray initial = getInitialStartingVector(vc);
        ListOfFloatArrays newNormVector = FinalScale.scaleToTargetVector(ba, matrixSize, initial, stem);
        if (newNormVector != null) {
            return ba.normalizeVectorByScaleFactor(newNormVector);
        } else {
            return null;
        }
    }
}