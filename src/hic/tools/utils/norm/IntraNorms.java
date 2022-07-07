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
import hic.tools.utils.bigarray.BigContactArrayCreator;
import hic.tools.utils.bigarray.BigContactList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.Map;
import java.util.Set;

public class IntraNorms {
    public static void getAllTheNorms(Dataset ds, HiCZoom zoom, int resolutionCutoffToSaveRAM,
                                      NormVectorsContainer container,
                                      boolean weShouldBuildVC, boolean weShouldBuildVCSqrt, boolean weShouldBuildScale,
                                      Map<NormalizationType, Integer> resolutionsToBuildTo,
                                      Set<Chromosome> scaleBPFailChroms) {

        for (Chromosome chrom : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chrom, chrom);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(zoom);
            if (zd == null) continue;

            if (HiCGlobals.printVerboseComments) {
                System.out.println("Now Doing " + chrom.getName());
            }

            BigContactList ba;
            if (zoom.getBinSize() < resolutionCutoffToSaveRAM) {
                ba = BigContactArrayCreator.createLocalVersionFromZD(zd);
            } else {
                ba = BigContactArrayCreator.createFromZD(zd);
            }
            matrix.clearCacheForZoom(zoom);

            NormalizationCalculations nc = new NormalizationCalculations(ba, zoom.getBinSize());

            if (weShouldBuildVC || weShouldBuildVCSqrt || weShouldBuildScale) {
                boolean saveVC = weShouldBuildVC && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.VC);
                boolean saveVCSqrt = weShouldBuildVCSqrt && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.VC_SQRT);
                boolean saveScale = weShouldBuildScale && zoom.getBinSize() >= resolutionsToBuildTo.get(NormalizationHandler.SCALE);

                if (saveVC || saveVCSqrt || saveScale) {
                    buildTheNorms(saveVC, saveVCSqrt, saveScale, chrom, nc, zoom, scaleBPFailChroms, container);
                }
            }

            nc = null;
            ba.clear();
        }
    }

    private static void buildTheNorms(boolean saveVC, boolean saveVCSqrt, boolean saveScale, Chromosome chrom,
                                      NormalizationCalculations nc, HiCZoom zoom,
                                      Set<Chromosome> scaleBPFailChroms, NormVectorsContainer container) {

        final int chrIdx = chrom.getIndex();
        ListOfFloatArrays vc = nc.computeVC();
        String stem = "NORM_" + chrIdx + "_" + zoom.getBinSize();

        if (saveScale) {
            if (!scaleBPFailChroms.contains(chrom)) {
                ListOfFloatArrays scale = nc.computeSCALE(vc, stem);
                if (scale == null) {
                    scaleBPFailChroms.add(chrom);
                } else {
                    nc.fixBySumFactor(scale);
                    container.add(chrom, new FloatNormVector(NormalizationHandler.SCALE, chrom.getIndex(), zoom, scale));
                    //updateExpectedValueCalculationForChr(chrIdx, nc, scale, NormalizationHandler.SCALE, zoom, ba, evSCALE, normVectorBuffers, normVectorIndices);
                }
            }
        }
        if (saveVCSqrt) {
            ListOfFloatArrays vcSqrt = new ListOfFloatArrays(vc.getLength());
            for (int i = 0; i < vc.getLength(); i++) {
                vcSqrt.set(i, (float) Math.sqrt(vc.get(i)));
            }
            nc.fixBySumFactor(vcSqrt);
            container.add(chrom, new FloatNormVector(NormalizationHandler.VC_SQRT, chrom.getIndex(), zoom, vcSqrt));
            //updateExpectedValueCalculationForChr(chrIdx, nc, vcSqrt, NormalizationHandler.VC_SQRT, zoom, ba, evVCSqrt, normVectorBuffers, normVectorIndices);
        }
        if (saveVC) {
            nc.fixBySumFactor(vc);
            container.add(chrom, new FloatNormVector(NormalizationHandler.VC, chrom.getIndex(), zoom, vc));
            //updateExpectedValueCalculationForChr(chrIdx, nc, vc, NormalizationHandler.VC, zoom, ba, evVC, normVectorBuffers, normVectorIndices);
        }
    }
}
