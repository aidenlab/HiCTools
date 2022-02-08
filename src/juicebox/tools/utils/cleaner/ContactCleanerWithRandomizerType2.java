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

package juicebox.tools.utils.cleaner;

import javastraw.reader.basics.ChromosomeHandler;
import juicebox.tools.utils.original.FragmentCalculation;

import java.util.List;
import java.util.Random;

public class ContactCleanerWithRandomizerType2 extends ContactCleanerWithRandomizer {
    private final List<FragmentCalculation> fragmentCalculations;

    public ContactCleanerWithRandomizerType2(ChromosomeHandler chromosomeHandler,
                                             List<FragmentCalculation> fragmentCalculationsForRandomization,
                                             Random random) {
        super(chromosomeHandler, null, random);
        this.fragmentCalculations = fragmentCalculationsForRandomization;
    }

    @Override
    protected void randomizePositions(int chr1, int chr2, int frag1, int frag2, int bp1, int bp2) {
        FragmentCalculation fragMap1 = findFragMap(fragmentCalculations, handler.getChromosomeFromIndex(chr1).getName(), bp1, frag1);
        FragmentCalculation fragMap2 = findFragMap(fragmentCalculations, handler.getChromosomeFromIndex(chr2).getName(), bp2, frag2);

        if (fragMap1 == null && fragMap2 == null) {
            //noMapFoundCount += 1;
            return;
        } else if (fragMap1 != null && fragMap2 != null && fragMap1 != fragMap2) {
            //mapDifferentCount += 1;
            return;
        }

        FragmentCalculation fragMapToUse;
        if (fragMap1 != null) {
            fragMapToUse = fragMap1;
        } else {
            fragMapToUse = fragMap2;
        }

        bp1 = randomizePos(fragMapToUse, handler.getChromosomeFromIndex(chr1).getName(), frag1);
        bp2 = randomizePos(fragMapToUse, handler.getChromosomeFromIndex(chr2).getName(), frag2);
    }
}
