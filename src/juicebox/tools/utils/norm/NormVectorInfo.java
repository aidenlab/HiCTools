/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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

package juicebox.tools.utils.norm;

import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.NormalizationType;
import org.broad.igv.tdf.BufferedByteWriter;

import java.util.List;
import java.util.Map;

class NormVectorInfo {

    private final List<BufferedByteWriter> normVectorBuffers;
    private final List<NormalizationVectorIndexEntry> normVectorIndices;
    private final Map<String, ExpectedValueFunction> expectedValueFunctionMap;
    private final Map<NormalizationType, Map<String, NormalizationVector>> normalizationVectorsMap;

    NormVectorInfo(Map<NormalizationType, Map<String, NormalizationVector>> normalizationVectorsMap, List<BufferedByteWriter> normVectorBuffers, List<NormalizationVectorIndexEntry> normVectorIndices,
                   Map<String, ExpectedValueFunction> expectedValueFunctionMap) {
        this.normalizationVectorsMap = normalizationVectorsMap;
        this.normVectorBuffers = normVectorBuffers;
        this.normVectorIndices = normVectorIndices;
        this.expectedValueFunctionMap = expectedValueFunctionMap;

    }

    public List<BufferedByteWriter> getNormVectorBuffers() {
        return normVectorBuffers;
    }

    public List<NormalizationVectorIndexEntry> getNormVectorIndices() {
        return normVectorIndices;
    }

    public Map<String, ExpectedValueFunction> getExpectedValueFunctionMap() {
        return expectedValueFunctionMap;
    }

    public Map<NormalizationType, Map<String, NormalizationVector>> getNormalizationVectorsMap() {
        return normalizationVectorsMap;
    }
}