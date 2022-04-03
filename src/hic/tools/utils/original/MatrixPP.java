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

package hic.tools.utils.original;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;

import java.io.File;
import java.io.IOException;
import java.util.Map;

public class MatrixPP {

    private static final int INTRA_CUTOFF = 500;
    private static final int INTER_CUTOFF = 5000;
    private final int MAX_SQRT = (int) Math.sqrt(Integer.MAX_VALUE);
    private final int chr1Idx;
    private final int chr2Idx;
    private final MatrixZoomDataPP[] zoomData;

    /**
     * Constructor for creating a matrix and initializing zoomed data at predefined resolution scales.  This
     * constructor is used when parsing alignment files.
     * c
     *
     * @param chr1Idx             Chromosome 1
     * @param chromosomeHandler
     * @param bpBinSizes
     * @param chr2Idx             Chromosome 2
     */
    public MatrixPP(int chr1Idx, int chr2Idx, ChromosomeHandler chromosomeHandler, int[] bpBinSizes,
                    int countThreshold, int v9DepthBase, int BLOCK_CAPACITY) {
        this.chr1Idx = chr1Idx;
        this.chr2Idx = chr2Idx;

        int nResolutions = bpBinSizes.length;

        zoomData = new MatrixZoomDataPP[nResolutions];

        int zoom = 0; //
        for (int idx = 0; idx < bpBinSizes.length; idx++) {
            int binSize = bpBinSizes[zoom];
			Chromosome chrom1 = chromosomeHandler.getChromosomeFromIndex(chr1Idx);
			Chromosome chrom2 = chromosomeHandler.getChromosomeFromIndex(chr2Idx);
	
			// Size block (submatrices) to be ~500 bins wide.

            long len = Math.max(chrom1.getLength(), chrom2.getLength());
            // for now, this will not be a long
            int nBins = (int) (len / binSize + 1);   // Size of chrom in bins
            int nColumns;
            if (chrom1.equals(chrom2)) {
                nColumns = getNumColumnsFromNumBins(nBins, binSize, INTRA_CUTOFF);
            } else {
                nColumns = getNumColumnsFromNumBins(nBins, binSize, INTER_CUTOFF);
            }
            zoomData[idx] = new MatrixZoomDataPP(chrom1, chrom2, binSize, nColumns, zoom, countThreshold, v9DepthBase, BLOCK_CAPACITY);
            zoom++;

        }
    }

    /**
     * Constructor for creating a matrix with a single zoom level at a specified bin size.  This is provided
     * primarily for constructing a whole-genome view.
     *
     * @param chr1Idx Chromosome 1
     * @param chr2Idx Chromosome 2
     * @param binSize Bin size
     */
    MatrixPP(int chr1Idx, int chr2Idx, int binSize, int blockColumnCount, ChromosomeHandler chromosomeHandler,
             int countThreshold, int v9DepthBase) {
        this.chr1Idx = chr1Idx;
        this.chr2Idx = chr2Idx;
        zoomData = new MatrixZoomDataPP[1];
        zoomData[0] = new MatrixZoomDataPP(chromosomeHandler.getChromosomeFromIndex(chr1Idx), chromosomeHandler.getChromosomeFromIndex(chr2Idx),
                binSize, blockColumnCount, 0, countThreshold, v9DepthBase);

    }

    private int getNumColumnsFromNumBins(int nBins, int binSize, int cutoff) {
        int nColumns = nBins / Preprocessor.BLOCK_SIZE + 1;
        if (binSize < cutoff) {
            long numerator = (long) nBins * binSize;
            long denominator = (long) Preprocessor.BLOCK_SIZE * cutoff;
            nColumns = (int) (numerator / denominator) + 1;
        }
        return Math.min(nColumns, MAX_SQRT - 1);
    }


    String getKey() {
        return "" + chr1Idx + "_" + chr2Idx;
    }


    public void incrementCount(int pos1, int pos2, int frag1, int frag2, float score, Map<String, ExpectedValueCalculation> expectedValueCalculations, File tmpDir) throws IOException {
        for (MatrixZoomDataPP aZoomData : zoomData) {
            aZoomData.incrementCount(pos1, pos2, score, expectedValueCalculations, tmpDir);
        }
    }

    public void parsingComplete() {
        for (MatrixZoomDataPP zd : zoomData) {
            if (zd != null) // fragment level could be null
                zd.parsingComplete();
        }
    }

    int getChr1Idx() {
        return chr1Idx;
    }

    int getChr2Idx() {
        return chr2Idx;
    }

    MatrixZoomDataPP[] getZoomData() {
        return zoomData;
    }

    /**
     * used by multithreaded code
     */
    public void mergeMatrices(MatrixPP otherMatrix) {
        if (otherMatrix != null) {
            for (MatrixZoomDataPP aZoomData : zoomData) {
                for (MatrixZoomDataPP bZoomData : otherMatrix.zoomData) {
                    if (aZoomData.getZoom() == bZoomData.getZoom()) {
                        aZoomData.mergeMatrices(bZoomData);
                    }
                }
            }
        }
    }
}
