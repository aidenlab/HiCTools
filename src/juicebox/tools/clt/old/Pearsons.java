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

import htsjdk.tribble.util.LittleEndianOutputStream;
import javastraw.matrices.BasicMatrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import juicebox.HiCGlobals;
import juicebox.tools.clt.CommandLineParser;
import juicebox.tools.clt.JuiceboxCLT;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Class for calculating Pearsons (separated out from Dump)
 * @author Neva Durand
 */
public class Pearsons extends JuiceboxCLT {

    private static final int BLOCK_TILE = 500;
    private String ofile = null;
    private HiCZoom.HiCUnit unit = null;
    private int binSize = 0;
    private Chromosome chromosome1;


    public Pearsons() {
        super(getBasicUsage() + "\n\t-p, --pearsons_all_resolutions: calculate Pearson's at all resolutions");
    }

    public static String getBasicUsage(){
        return "pearsons [-p] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]";
    }

    @Override
    public void readArguments(String[] args, CommandLineParser parser) {
        if (args.length != 7 && args.length != 6) {
            printUsageAndExit();
        }

        //HiCGlobals.MAX_PEARSON_ZOOM = 500000;
        setDatasetAndNorm(args[2], args[1], true);
        ChromosomeHandler chromosomeHandler = dataset.getChromosomeHandler();

        if (chromosomeHandler.doesNotContainChromosome(args[3])) {
            System.err.println("Unknown chromosome: " + args[3]);
            System.exit(18);
        }
        chromosome1 = chromosomeHandler.getChromosomeFromName(args[3]);

        try {
            unit = HiCZoom.valueOfUnit(args[4]);
        } catch (IllegalArgumentException error) {
            System.err.println("Unit must be in BP or FRAG.");
            System.exit(20);
        }

        String binSizeSt = args[5];

        try {
            binSize = Integer.parseInt(binSizeSt);
        } catch (NumberFormatException e) {
            System.err.println("Integer expected for bin size.  Found: " + binSizeSt + ".");
            System.exit(21);
        }

        if ((unit == HiCZoom.HiCUnit.BP && binSize < HiCGlobals.MAX_PEARSON_ZOOM) ||
                (unit == HiCZoom.HiCUnit.FRAG && binSize < HiCGlobals.MAX_PEARSON_ZOOM / 1000)) {
            System.out.println("Pearson's and Eigenvector are not calculated for high resolution datasets");
            System.out.println("To override this limitation, send in the \"-p\" flag.");
            System.exit(0);
        }

        if (args.length == 7) {
            ofile = args[6];
        }

    }

    @Override
    public void run() {
        HiCZoom zoom = new HiCZoom(unit, binSize);
        MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chromosome1, chromosome1, zoom);
        if (zd == null) {
            System.err.println("No reads in " + chromosome1);
            System.err.println("Unknown resolution: " + zoom);
            System.err.println("This data set has the following bin sizes (in bp): ");
            for (int zoomIdx = 0; zoomIdx < dataset.getNumberZooms(HiCZoom.HiCUnit.BP); zoomIdx++) {
                System.err.print(dataset.getZoom(HiCZoom.HiCUnit.BP, zoomIdx).getBinSize() + " ");
            }
            System.err.println("\nand the following bin sizes (in frag): ");
            for (int zoomIdx = 0; zoomIdx < dataset.getNumberZooms(HiCZoom.HiCUnit.FRAG); zoomIdx++) {
                System.err.print(dataset.getZoom(HiCZoom.HiCUnit.FRAG, zoomIdx).getBinSize() + " ");
            }
            System.exit(13);
        }
        ExpectedValueFunction df = dataset.getExpectedValuesOrExit(zd.getZoom(), norm, chromosome1, true, true);

        BasicMatrix pearsons = zd.getPearsons(df);
        if (pearsons == null) {
            System.err.println("Pearson's not available at zoom " + zoom  + ". For high resolution, try again with -p");
            System.exit(15);
        }

        LittleEndianOutputStream les = null;
        PrintWriter txtWriter = null;
        if (ofile != null) {
            try {
                if (ofile.endsWith(".bin")) {
                    BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(ofile));
                    les = new LittleEndianOutputStream(bos);
                } else {
                    txtWriter = new PrintWriter(new FileOutputStream(ofile));
                }
            } catch (IOException error) {
                System.err.println("Cannot write to " + ofile);
                error.printStackTrace();
                System.exit(22);
            }
        }
        else {
            txtWriter = new PrintWriter(System.out);
        }

        if (les == null) {
            int dim = pearsons.getRowDimension();
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    float output = pearsons.getEntry(i, j);
                    txtWriter.print(output + " ");
                }
                txtWriter.println();
            }
            txtWriter.flush();
        }
        else {
            try {
                int dim = pearsons.getRowDimension();
                writeHeader(les, dim, pearsons.getLowerValue(), pearsons.getUpperValue());
                int block_side = (int) Math.ceil((float) dim / (float) BLOCK_TILE);
                for (int i = 0; i < block_side; i++) {
                    int block_row_start = i * BLOCK_TILE;
                    int block_row_end = Math.min(block_row_start + BLOCK_TILE, dim);
                    int row_len = block_row_end - block_row_start;
                    for (int j = 0; j < block_side; j++) {
                        int block_col_start = j * BLOCK_TILE;
                        int block_col_end = Math.min(block_col_start + BLOCK_TILE, dim);
                        int col_len = block_col_end - block_col_start;
                        for (int ui = 0; ui < row_len; ui++) {
                            for (int uj = 0; uj < col_len; uj++) {
                                int now_i = ui + block_row_start;
                                int now_j = uj + block_col_start;
                                float output = pearsons.getEntry(now_i, now_j);
                                les.writeFloat(output);
                            }
                        }

                    }
                }
                les.close();
            }
            catch (IOException error) {
                System.err.println("Problem when writing Pearson's");
                error.printStackTrace();
                System.exit(1);
            }
        }
    }

    private void writeHeader(LittleEndianOutputStream les, int dim, float lower, float upper) throws IOException {

        // Magic number - 4 bytes
        les.writeByte('h');
        les.writeByte('i');
        les.writeByte('c');
        les.writeByte(0);

        // Version number
        les.writeInt(1);

        // Genome --
        les.writeString(dataset.getGenomeId());

        // Chromosomes
        les.writeString(chromosome1.getName());
        les.writeString(chromosome1.getName());

        // Resolution (bin size)
        les.writeInt(binSize);

        // Statistics, other attributes
        les.writeFloat(lower);  // this is supposed to be lower quartile
        les.writeFloat(upper);  // this is supposed to be upper quartile
        les.writeInt(dim);  // # rows
        les.writeInt(dim);  // # cols
        les.writeInt(BLOCK_TILE);
    }
}
