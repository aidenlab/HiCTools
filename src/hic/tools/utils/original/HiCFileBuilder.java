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

import hic.tools.utils.mnditerator.ReadPairFilter;
import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.*;
import java.util.zip.Deflater;

abstract public class HiCFileBuilder {

    protected static final int VERSION = 9;
    protected static final int BLOCK_SIZE = 1000;
    public static int BLOCK_CAPACITY = 1000;
    protected final ChromosomeHandler chromosomeHandler;
    protected final File outputFile;
    protected final Map<String, IndexEntry> matrixPositions = new LinkedHashMap<>();
    protected final Deflater compressor = WriterUtils.getDefaultCompressor();
    protected int v9DepthBase = 2;
    protected Map<String, Integer> chromosomeIndexes = new Hashtable<>();
    protected String genomeId;
    protected final LittleEndianOutputStream[] losArray = new LittleEndianOutputStream[1];
    protected long masterIndexPosition;
    protected int countThreshold = 0;
    protected int mapqThreshold = 0;
    protected boolean diagonalsOnly = false;
    protected String statsFileName = null;
    protected String graphFileName = null;
    protected Set<String> includedChromosomes;
    protected ReadPairFilter filter = null;
    protected int[] bpBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000};
    protected int numResolutions;
    protected double hicFileScalingFactor = 1;
    protected File tmpDir = UNIXTools.makeDir(new File("temp_folder"));

    public HiCFileBuilder(File outputFile, String genomeId, ChromosomeHandler chromosomeHandler, double hicFileScalingFactor) {
        this.genomeId = genomeId;
        this.outputFile = outputFile;
        this.chromosomeHandler = chromosomeHandler;
        for (int i = 0; i < chromosomeHandler.size(); i++) {
            chromosomeIndexes.put(chromosomeHandler.getChromosomeFromIndex(i).getName(), i);
        }
        if (hicFileScalingFactor > 0) {
            this.hicFileScalingFactor = hicFileScalingFactor;
        }
    }

    public void setCountThreshold(int countThreshold) {
        this.countThreshold = countThreshold;
    }

    public void setV9DepthBase(int v9DepthBase) {
        if (v9DepthBase > 1 || v9DepthBase < 0) {
            this.v9DepthBase = v9DepthBase;
        }
    }

    public void setMapqThreshold(int mapqThreshold) {
        this.mapqThreshold = mapqThreshold;
    }

    public void setDiagonalsOnly(boolean diagonalsOnly) {
        this.diagonalsOnly = diagonalsOnly;
    }

    public void setIncludedChromosomes(Set<String> includedChromosomes) {
        if (includedChromosomes != null && includedChromosomes.size() > 0) {
            this.includedChromosomes = Collections.synchronizedSet(new HashSet<>());
            for (String name : includedChromosomes) {
                this.includedChromosomes.add(chromosomeHandler.cleanUpName(name));
            }
        }
    }

    public void setGraphFile(String graphFileName) {
        this.graphFileName = graphFileName;
    }

    public void setGenome(String genome) {
        if (genome != null) {
            this.genomeId = genome;
        }
    }

    public void setResolutions(List<String> resolutions) {
        if (resolutions != null) {
            ArrayList<Integer> bpResolutions = new ArrayList<>();
            for (String str : resolutions) {
                try {
                    int myInt = Integer.parseInt(str);
                    bpResolutions.add(myInt);
                } catch (NumberFormatException exception) {
                    System.err.println("Resolution improperly formatted. It must be in the form of a number, such as 1000000 for 1 million BP");
                    System.err.println("Fragment resolutions are deprecated and require using an older jar.");
                    System.exit(1);
                }
            }

            boolean resolutionsSet = false;
            if (bpResolutions.size() > 0) {
                resolutionsSet = true;
                Collections.sort(bpResolutions);
                Collections.reverse(bpResolutions);
                int[] bps = new int[bpResolutions.size()];
                for (int i = 0; i < bps.length; i++) {
                    bps[i] = bpResolutions.get(i);
                }
                bpBinSizes = bps;
            } else {
                bpBinSizes = new int[0];
            }
            if (!resolutionsSet) {
                System.err.println("No valid resolutions sent in");
                System.exit(1);
            }
        }
    }

    public void setFilter(ReadPairFilter.Type type) {
        if (type != null) {
            this.filter = new ReadPairFilter(type);
        }
    }


    public void setTmpdir(String tmpDirName) {
        if (tmpDirName != null) {
            this.tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                System.err.println("Tmp directory does not exist: " + tmpDirName);
                if (outputFile != null) outputFile.deleteOnExit();
                System.exit(59);
            }
        }
    }

    public void setStatisticsFile(String statsOption) {
        statsFileName = statsOption;
    }
}
