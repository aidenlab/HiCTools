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

import hic.tools.utils.largelists.BigFloatsArray;
import hic.tools.utils.largelists.BigShortsArray;
import hic.tools.utils.original.ExpectedValueCalculation;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;

public interface BigContactList {
    void updateGenomeWideExpected(int chrIdx, ListOfFloatArrays vector, ExpectedValueCalculation exp);

    long getMatrixSize();

    void clear();

    ListOfFloatArrays getRowSums();

    double[] getNormMatrixSumFactor(ListOfFloatArrays norm);

    ListOfFloatArrays normalizeVectorByScaleFactor(ListOfFloatArrays newNormVector);

    ListOfIntArrays getNumNonZeroInRows();

    BigFloatsArray parSparseMultiplyAcrossLists(BigShortsArray one, long matrixSize);

    BigFloatsArray parSparseMultiplyAcrossLists(BigFloatsArray one, long matrixSize);

    void clearIntraAndShiftInter();
}
