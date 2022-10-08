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

package hic.tools.clt.old;

import hic.HiCGlobals;
import hic.tools.clt.CommandLineParser;
import hic.tools.clt.JuiceboxCLT;
import javastraw.reader.Dataset;
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;


public class DumpNorm extends JuiceboxCLT {
    private String file, output;

    public DumpNorm() {
        super(getBasicUsage());
    }

    public static String getBasicUsage() {
        return "dump-norms <input.hic> <output_>";
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {
        if (parser.getHelpOption() || args.length != 3) {
            printUsageAndExit();
        }

        HiCGlobals.normThreads = updateNumberOfCPUThreads(parser, 10);

        file = args[1];
        output = args[2];
    }

    @Override
    public void run() {
        HiCGlobals.allowDynamicBlockIndex = false;

        try {
            DatasetReaderV2 reader = new DatasetReaderV2(file, false, false);
            Dataset ds = reader.read();
            HiCGlobals.verifySupportedHiCFileVersion(reader.getVersion());
            for (Chromosome chrom : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
                for (HiCZoom zoom : ds.getBpZooms()) {
                    for (NormalizationType norm : ds.getNormalizationHandler().getDefaultSetForHiCFileBuilding()) {
                        String outString = output + "_" + chrom.getName() + "_" + zoom.getBinSize() + "_" + norm.getLabel() + ".txt";
                        TestNorm.export(ds.getNormalizationVector(chrom.getIndex(), zoom, norm).getData(), outString);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}