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

import javastraw.tools.ParallelizationTools;
import juicebox.HiCGlobals;

import java.util.concurrent.atomic.AtomicInteger;

public class ParallelVectorTools {
    public static void process(long matrixSize, VectorFunction function) {

        long[] cutoffs = getCutoffs(matrixSize);
        int numThreads = HiCGlobals.numCPUMatrixThreads;

        AtomicInteger index = new AtomicInteger();
        ParallelizationTools.launchParallelizedCode(numThreads, () -> {

            int i = index.getAndIncrement();
            while (i < numThreads) {
                for (long p = cutoffs[i]; p < cutoffs[i + 1]; p++) {
                    function.use(p);
                }
                i = index.getAndIncrement();
            }
        });
    }

    private static long[] getCutoffs(long matrixSize) {
        int n = HiCGlobals.numCPUMatrixThreads;
        long[] bounds = new long[n + 1];
        for (int z = 0; z < n; z++) {
            bounds[z] = (matrixSize * z / n);
        }
        bounds[0] = 0;
        bounds[n] = matrixSize;
        return bounds;
    }
}
