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
import juicebox.tools.utils.mnditerator.AlignmentPair;
import juicebox.tools.utils.original.FragmentCalculation;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class ContactCleanerWithRandomizer extends ContactCleaner {
    private final FragmentCalculation fragmentCalculation;
    private final Random random;

    public ContactCleanerWithRandomizer(ChromosomeHandler chromosomeHandler, FragmentCalculation fragmentCalculation,
                                        Random random) {
        super(chromosomeHandler);
        this.fragmentCalculation = fragmentCalculation;
        this.random = random;
    }

    public static FragmentCalculation findFragMap(List<FragmentCalculation> maps, String chr, int bp, int frag) {
        //potential maps that this strand could come from
        ArrayList<FragmentCalculation> mapsFound = new ArrayList<>();
        for (FragmentCalculation fragmentCalculation : maps) {
            int low = 1;
            int high = 1;

            if (frag > fragmentCalculation.getNumberFragments(chr)) {
                // definitely not this restriction site file for certain
                continue;
            }

            try {
                if (frag == 0) {
                    high = fragmentCalculation.getSites(chr)[frag];
                } else if (frag == fragmentCalculation.getNumberFragments(chr)) {
                    high = fragmentCalculation.getSites(chr)[frag - 1];
                    low = fragmentCalculation.getSites(chr)[frag - 2];
                } else {
                    high = fragmentCalculation.getSites(chr)[frag];
                    low = fragmentCalculation.getSites(chr)[frag - 1];
                }
            } catch (Exception e) {
                e.printStackTrace();
                System.out.printf("fragment: %d, number of frags: %d%n", frag, fragmentCalculation.getNumberFragments(chr));

            }

            // does bp fit in this range?
            if (bp >= low && bp <= high) {
                mapsFound.add(fragmentCalculation);
            }
        }
        if (mapsFound.size() == 1) {
            return mapsFound.get(0);
        }
        return null;
    }

    @Override
    public void updateLatestContact(AlignmentPair pair) {
        super.updateLatestContact(pair);
        randomizePositions(chr1, chr2, frag1, frag2, bp1, bp2);
    }

    protected void randomizePositions(int chr1, int chr2, int frag1, int frag2, int bp1, int bp2) {
        FragmentCalculation fragMapToUse = fragmentCalculation;
        bp1 = randomizePos(fragMapToUse, handler.getChromosomeFromIndex(chr1).getName(), frag1);
        bp2 = randomizePos(fragMapToUse, handler.getChromosomeFromIndex(chr2).getName(), frag2);
    }

    protected int randomizePos(FragmentCalculation fragmentCalculation, String chr, int frag) {
        int low = 1;
        int high = 1;
        if (frag == 0) {
            high = fragmentCalculation.getSites(chr)[frag];
        } else if (frag >= fragmentCalculation.getNumberFragments(chr)) {
            high = fragmentCalculation.getSites(chr)[frag - 1];
            low = fragmentCalculation.getSites(chr)[frag - 2];
        } else {
            high = fragmentCalculation.getSites(chr)[frag];
            low = fragmentCalculation.getSites(chr)[frag - 1];
        }
        return random.nextInt(high - low + 1) + low;
    }
}
