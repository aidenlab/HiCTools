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

package juicebox.tools.utils.norm.scale;

import javastraw.reader.datastructures.ListOfFloatArrays;
import juicebox.tools.utils.bigarray.BigArray;

public class ScaleHandler {
    public static ListOfFloatArrays scale(BigArray ic, int resolution, long matrixSize,
                                          ListOfFloatArrays initialGuess) {
        /*
        BigArray ic2 = ic;
        if (DistanceFilteredIteratorContainer.getUseFilterDistance()) {
            ic2 = new DistanceFilteredIteratorContainer(ic, resolution);
        }
        */
        return FinalScale.scaleToTargetVector(ic, matrixSize, initialGuess);
    }


    public static ListOfFloatArrays mmbaScaleToVector(BigArray ic, int resolution, long matrixSize,
                                                      ListOfFloatArrays initialGuess) {
        ListOfFloatArrays newNormVector = scale(ic, resolution, matrixSize, initialGuess);
        return ic.normalizeVectorByScaleFactor(newNormVector);
    }

}
