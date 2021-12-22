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

package juicebox.tools.clt.old;

import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import juicebox.tools.clt.CommandLineParser;
import juicebox.tools.clt.JuiceboxCLT;
import juicebox.tools.utils.original.mnditerator.AsciiToBinConverter;

public class PairsToBin extends JuiceboxCLT {

    private String ifile, ofile, genomeId;

    public PairsToBin() {
        super("pairsToBin <input_mnd> <output_mnd_binary> <genomeID>");
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }
        ifile = args[1];
        ofile = args[2];
        genomeId = args[3];
    }

    @Override
    public void run() {
        ChromosomeHandler chromosomeHandler = ChromosomeTools.loadChromosomes(genomeId);
        try {
            AsciiToBinConverter.convert(ifile, ofile, chromosomeHandler);
        } catch (Exception e) {
            System.err.println("Unable to convert from ascii to bin");
            e.printStackTrace();
        }
    }
}
