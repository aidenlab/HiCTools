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

import javastraw.reader.basics.Chromosome;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class NormVectorsContainer {
    private final Map<NormalizationType, Map<Chromosome, FloatNormVector>> allData = new HashMap<>();
    private final List<NormalizationType> gwNormalizations;
    private final List<NormalizationType> gwInterNormalizations;

    public NormVectorsContainer(List<NormalizationType> normalizationsToBuild, Map<NormalizationType, Integer> resolutionsToBuildTo, HiCZoom zoom) {
        gwNormalizations = GWNorms.getGWNorms(normalizationsToBuild, resolutionsToBuildTo, zoom);
        gwInterNormalizations = GWNorms.getInterNorms(normalizationsToBuild, resolutionsToBuildTo, zoom);
    }

    public static List<NormalizationType> sortedNorms() {
        List<NormalizationType> norms0 = new ArrayList<>();
        norms0.add(NormalizationHandler.GW_SCALE);
        norms0.add(NormalizationHandler.GW_VC);
        norms0.sort(Comparator.comparing(NormalizationType::getLabel));

        List<NormalizationType> norms1 = new ArrayList<>();
        norms1.add(NormalizationHandler.INTER_SCALE);
        norms1.add(NormalizationHandler.INTER_VC);
        norms1.sort(Comparator.comparing(NormalizationType::getLabel));
        norms0.addAll(norms1);

        norms0.add(NormalizationHandler.VC);
        norms0.add(NormalizationHandler.VC_SQRT);
        norms0.add(NormalizationHandler.SCALE);

        return norms0;
    }

    public boolean hasNoGenomewideNorms() {
        return gwNormalizations.isEmpty() && gwInterNormalizations.isEmpty();
    }

    public boolean useGWIntra() {
        return !gwNormalizations.isEmpty();
    }

    public List<NormalizationType> getGenomewideNorms() {
        return gwNormalizations;
    }

    public List<NormalizationType> getGenomewideInterNorms() {
        return gwInterNormalizations;
    }

    public synchronized void put(NormalizationType normType, Map<Chromosome, FloatNormVector> wgVectors) {
        allData.put(normType, wgVectors);
    }

    public synchronized void add(Chromosome chrom, FloatNormVector vector) {
        NormalizationType norm = vector.getNormType();
        if (!allData.containsKey(norm)) {
            allData.put(norm, new HashMap<>());
        }
        allData.get(norm).put(chrom, vector);
    }

    public synchronized Map<Chromosome, FloatNormVector> get(NormalizationType norm) {
        return allData.get(norm);
    }

    public synchronized boolean containsNorm(NormalizationType norm) {
        return allData.containsKey(norm);
    }

    public synchronized void clear() {
        for (NormalizationType norm : allData.keySet()) {
            allData.get(norm).clear();
        }
        allData.clear();
    }

    public synchronized Set<NormalizationType> getNorms() {
        return allData.keySet();
    }
}
