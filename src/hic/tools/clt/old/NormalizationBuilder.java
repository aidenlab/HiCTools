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

package hic.tools.clt.old;

import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

public class NormalizationBuilder {

    public static Integer getIdealResolutionLimit(NormalizationType normalizationType) {
        if (NormalizationHandler.isGenomeWideNorm(normalizationType)) {
            return 1000;
        } else {
            return 0;
        }
    }

    public static boolean usesVC(NormalizationType option) {
        return option.getLabel().contains("VC");
    }

    public static boolean usesSCALE(NormalizationType option) {
        return option.getLabel().contains("SCALE");
    }

    public static boolean isNONE(NormalizationType option) {
        return option.getLabel().equals("NONE");
    }
}
