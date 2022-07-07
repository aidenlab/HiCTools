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

import hic.HiCGlobals;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ParallelizationTools;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class NormalizationTools {
    public static Map<Chromosome, FloatNormVector> parCreateNormVectorMap(ChromosomeHandler chromosomeHandler,
                                                                          int resolution, ListOfFloatArrays vector,
                                                                          NormalizationType norm, HiCZoom zoom) {
        final Map<Chromosome, FloatNormVector> normVectorMap = new LinkedHashMap<>();

        final AtomicInteger index = new AtomicInteger(0);
        Chromosome[] chromosomes = chromosomeHandler.getChromosomeArrayWithoutAllByAll();
        long[] offsets = createOffsets(chromosomes, resolution);
        ParallelizationTools.launchParallelizedCode(HiCGlobals.normThreads, () -> {
            int i = index.getAndIncrement();
            while (i < (chromosomes).length) {
                Chromosome c1 = chromosomes[i];
                long offset = offsets[i];
                long chrBinned = c1.getLength() / resolution + 1;
                ListOfFloatArrays chrNV = new ListOfFloatArrays(chrBinned);
                for (long k = 0; k < chrNV.getLength(); k++) { // todo optimize a version with system.arraycopy and long
                    chrNV.set(k, vector.get(offset + k));
                }
                synchronized (normVectorMap) {
                    normVectorMap.put(c1, new FloatNormVector(norm, c1.getIndex(), zoom, chrNV));
                }
                i = index.getAndIncrement();
            }
        });

        return normVectorMap;
    }

    private static long[] createOffsets(Chromosome[] chromosomes, int resolution) {
        long[] offsets = new long[chromosomes.length];
        offsets[0] = 0L;
        for (int i = 0; i < chromosomes.length - 1; i++) {
            offsets[i + 1] = offsets[i] + (chromosomes[i].getLength() / resolution) + 1;
        }
        return offsets;
    }

    public static boolean checkIfInterDataAvailable(HiCZoom zoom, Dataset ds) {

        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i + 1; j < chromosomes.length; j++) {
                try {
                    Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[j]);
                    if (matrix == null) return false;
                    MatrixZoomData zd = matrix.getZoomData(zoom);
                    if (zd == null) return false;
                } catch (Exception e) {
                    return false;
                }
            }
        }
        return true;
    }
}
