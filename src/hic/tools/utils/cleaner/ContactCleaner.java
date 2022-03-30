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

package hic.tools.utils.cleaner;

import hic.tools.utils.mnditerator.AlignmentPair;
import hic.tools.utils.original.ExpectedValueCalculation;
import hic.tools.utils.original.MatrixPP;
import javastraw.reader.basics.ChromosomeHandler;

import java.io.File;
import java.io.IOException;
import java.util.Map;

public class ContactCleaner {
    protected final ChromosomeHandler handler;
    protected int chr1, chr2, bp1, bp2, frag1, frag2;
    private float score;

    public ContactCleaner(ChromosomeHandler chromosomeHandler) {
        this.handler = chromosomeHandler;
    }

    public static int getWholeGenomePosition(int chr, int pos, ChromosomeHandler handler) {
        long len = 0;
        for (int i = 1; i < chr; i++) {
            len += handler.getChromosomeFromIndex(i).getLength();
        }
        len += pos;

        return (int) (len / 1000);
    }

    public void updateLatestContact(AlignmentPair pair) {
        if (isUpperTriangular(pair)) {
            bp1 = pair.getPos1();
            bp2 = pair.getPos2();
            frag1 = pair.getFrag1();
            frag2 = pair.getFrag2();
            chr1 = pair.getChr1();
            chr2 = pair.getChr2();
        } else {
            bp1 = pair.getPos2();
            bp2 = pair.getPos1();
            frag1 = pair.getFrag2();
            frag2 = pair.getFrag1();
            chr1 = pair.getChr2();
            chr2 = pair.getChr1();
        }
        score = pair.getScore();

        bp1 = ensureFitInChromosomeBounds(bp1, chr1);
        bp2 = ensureFitInChromosomeBounds(bp2, chr2);
    }

    private boolean isUpperTriangular(AlignmentPair pair) {
        boolean isIntra = pair.getChr1() == pair.getChr2();
        boolean isUpperTriangular = pair.getPos1() <= pair.getPos2();
        boolean isInterUpperTriangular = pair.getChr1() < pair.getChr2();
        return (isIntra && isUpperTriangular) || isInterUpperTriangular;
    }

    protected int ensureFitInChromosomeBounds(int bp, int chrom) {
        if (bp < 0) {
            return 0;
        }
        long maxLength = handler.getChromosomeFromIndex(chrom).getLength();
        if (bp > maxLength) {
            return (int) maxLength;
        }
        return bp;
    }

    public void incrementCount(MatrixPP currentMatrix, Map<String, ExpectedValueCalculation> expectedValueCalculations,
                               File tmpDir) throws IOException {
        if (currentMatrix != null) {
            currentMatrix.incrementCount(bp1, bp2, frag1, frag2, score, expectedValueCalculations, tmpDir);
        }
    }

    public boolean doesntMatchCurrentBlock(int currentChr1, int currentChr2) {
        return !(currentChr1 == chr1 && currentChr2 == chr2);
    }

    public int getChr1() {
        return chr1;
    }

    public int getChr2() {
        return chr2;
    }

    public void incrementGWCount(MatrixPP wholeGenomeMatrix, Map<String, ExpectedValueCalculation> localExpectedValueCalculations, File tmpDir) throws IOException {
        int pos1 = getWholeGenomePosition(chr1, bp1, handler);
        int pos2 = getWholeGenomePosition(chr2, bp2, handler);
        wholeGenomeMatrix.incrementCount(pos1, pos2, pos1, pos2, score, localExpectedValueCalculations, tmpDir);
    }

    public String getCurrentMatrixName() {
        return handler.getChromosomeFromIndex(chr1).getName()
                + "-" + handler.getChromosomeFromIndex(chr2).getName();
    }
}
