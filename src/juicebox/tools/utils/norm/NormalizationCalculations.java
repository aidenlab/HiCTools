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

package juicebox.tools.utils.norm;

import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.type.NormalizationType;
import juicebox.tools.clt.old.NormalizationBuilder;
import juicebox.tools.utils.bigarray.BigContactArray;
import juicebox.tools.utils.largelists.BigFloatsArray;
import juicebox.tools.utils.norm.scale.ScaleHandler;

public class NormalizationCalculations {

    private final long matrixSize; // x and y symmetric
    //private boolean isEnoughMemory = false;
    private final BigContactArray ba;
    private final int resolution;

    public NormalizationCalculations(BigContactArray ba, int resolution) {
        this.ba = ba;
        this.matrixSize = ba.getMatrixSize();
        this.resolution = resolution;
    }

    public ListOfFloatArrays getNorm(NormalizationType normOption) {
        ListOfFloatArrays norm;
        if (NormalizationBuilder.usesVC(normOption)) {
            norm = computeVC();
        } else if (NormalizationBuilder.usesSCALE(normOption)) {
            norm = computeSCALE(computeVC());
        } else if (NormalizationBuilder.isNONE(normOption)) {
            return new ListOfFloatArrays(matrixSize, 1);
        } else {
            System.err.println("Not supported for normalization " + normOption);
            return null;
        }

        if (norm != null && norm.getLength() > 0) {
            double factor = getSumFactor(norm);
            norm.multiplyEverythingBy(factor);
        }
        return norm;
    }
    
    /**
     * Compute vanilla coverage norm, just the sum of the rows
     *
     * @return Normalization vector
     */
    ListOfFloatArrays computeVC() {
        return ba.getRowSums();
    }

    /**
     * Get the sum of the normalized matrix
     *
     * @param norm Normalization vector
     * @return Square root of ratio of original to normalized vector
     */
    public double getSumFactor(ListOfFloatArrays norm) {
        double[] normMatrixSums = getNormMatrixSumFactor(norm);
        return Math.sqrt(normMatrixSums[0] / normMatrixSums[1]);
    }
    
    public double[] getNormMatrixSumFactor(ListOfFloatArrays norm) {
        return ba.getNormMatrixSumFactor(norm);
    }

    private BigFloatsArray getInitialStartingVector(ListOfFloatArrays vc) {
        BigFloatsArray initial = new BigFloatsArray(vc.getLength());
        for (long i = 0; i < vc.getLength(); i++) {
            initial.set(i, (float) Math.sqrt(vc.get(i)));
        }
        return initial;
    }

    public ListOfFloatArrays computeSCALE(ListOfFloatArrays vc) {
        BigFloatsArray initial = getInitialStartingVector(vc);
        return ScaleHandler.mmbaScaleToVector(ba, matrixSize, initial);
    }

    /*public BigContactRecordList booleanBalancing() {
        ListOfFloatArrays rowsums = new ListOfFloatArrays(totSize, 0);
        Map<Float,Map<Long,Integer>> rowsumIndices = new HashMap<>();
        Map<Long, List<LinkedContactRecord>> rows = new HashMap<>();
        //Map<Long, RandomizedCollection> remainingContacts = new HashMap<>();
        List<Double> thresholds = new ArrayList<>(Arrays.asList(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,3.0,4.0,5.0));

        for (List<ContactRecord> localList : contactRecords) {
            for (ContactRecord cr : localList) {
                int x = cr.getBinX();
                int y = cr.getBinY();
                float value = cr.getCounts();
                rowsums.addTo(x, value);
                if (x != y) {
                    rowsums.addTo(y, value);

                }
                List<LinkedContactRecord> row1 = rows.get((long) x);
                List<LinkedContactRecord> row2 = rows.get((long) y);
                if (row1 == null) {
                    row1 = new ArrayList<>();
                    rows.put((long) x, row1);
                    //remainingContacts.put((long) x, new RandomizedCollection());
                }
                if (row2 == null) {
                    row2 = new ArrayList<>();
                    rows.put((long) y, row2);
                    //remainingContacts.put((long) y, new RandomizedCollection());
                }
                int xCurrentSize = row1.size();
                int yCurrentSize = row2.size();
                for (int i = 0; i < value; i++) {
                    row1.add(new LinkedContactRecord(cr, yCurrentSize+i));
                    //remainingContacts.get((long) x).insert(xCurrentSize);
                    if (x != y) {
                        row2.add(new LinkedContactRecord(cr, xCurrentSize+i));
                        //remainingContacts.get((long) y).insert(yCurrentSize);
                    }
                }


            }
        }
        System.out.println("loaded contacts for matrix: " + chr1 + "-" + chr2);

        int rowSumThreshold = (int) (1.0d / getSumFactor(rowsums));
        List<Float> sortedRowSums = new ArrayList<>();
        for (long i = 0; i < rowsums.getLength(); i++) {
            //Instant E = Instant.now();
            Map<Long,Integer> sumIndexList = rowsumIndices.get(rowsums.get(i));
            if (sumIndexList == null) {
                rowsumIndices.put(rowsums.get(i), new HashMap<>());
                sumIndexList = rowsumIndices.get(rowsums.get(i));
            }
            sumIndexList.put(i,1);
            sortedRowSums.add(rowsums.get(i));
            //Instant F = Instant.now();
            //System.err.println(Duration.between(E,F).toNanos());
        }

        System.out.println("mapped row sums to rows indices for matrix: " + chr1 + "-" + chr2);


        Collections.sort(sortedRowSums);

        Map<Float,Integer> SumMap = new HashMap<>();
        List<Float> SumKeys = new ArrayList<>(rowsumIndices.keySet());
        Collections.sort(SumKeys);
        int keyCounter = 0;
        for (float key : SumKeys) {
            //Collections.sort(rowsumIndices.get(key));
            while (sortedRowSums.get(keyCounter)!=key) {
                keyCounter++;
            }
            SumMap.put(key,keyCounter);
        }

        System.out.println("sorted row sums for matrix: " + chr1 + "-" + chr2);

        Instant A = Instant.now();
        Map<Long,Integer> currentRows = rowsumIndices.get(sortedRowSums.get(sortedRowSums.size()-1));
        long currentRow = currentRows.keySet().iterator().next();
        //System.out.println(currentRows.keySet().size() + " " + sortedRowSums.get(sortedRowSums.size()-1) + " " + currentRow + " " + rowsums.get(currentRow));
        Instant B = Instant.now();
        //System.err.println(Duration.between(A,B).toNanos());
        //System.err.println(rowSumThreshold + " " + getSumFactor(rowsums) + " " + rowsums.get(rowsums.getMaxRow()));
        Random randomNumberGenerator = new Random(0);
        while (rowsums.get(currentRow) > rowSumThreshold) {
            float currentRowSum = rowsums.get(currentRow);
            //Instant C = Instant.now();
            //List<Integer> removedContactList = new ArrayList<>(removedContacts.get(currentRow));
            //Collections.sort(removedContactList);
            //int randomContact = getRandomWithExclusion(randomNumberGenerator, rows.get(currentRow).size(), removedContactList);
            //removedContacts.get(currentRow).add(randomContact);
            //int randomContact = remainingContacts.get(currentRow).getRandom();
            int randomContact = (int) (rows.get(currentRow).size() * randomNumberGenerator.nextDouble());
            //remainingContacts.get(currentRow).remove(randomContact);
            long firstRow, secondRow;
            if (rows.get(currentRow).get(randomContact).getContactRecord().getBinX() == currentRow) {
                firstRow = (long) rows.get(currentRow).get(randomContact).getContactRecord().getBinX();
                secondRow = (long) rows.get(currentRow).get(randomContact).getContactRecord().getBinY();
            } else {
                firstRow = (long) rows.get(currentRow).get(randomContact).getContactRecord().getBinY();
                secondRow = (long) rows.get(currentRow).get(randomContact).getContactRecord().getBinX();
            }
            //System.out.println(currentRow + " " + firstRow + " " + secondRow);
            float firstRowSum = rowsums.get(firstRow);
            int symmetricRandomContactPosition = rows.get(firstRow).get(randomContact).getLink();
            long lastRow = rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinX() == firstRow? rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinY() : rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinX();
            int lastLink = rows.get(firstRow).get(rows.get(firstRow).size()-1).getLink();
            //System.out.println("initial values " + currentRow + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinX() + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinY() + " " + rows.get(firstRow).get(randomContact).getLink() +  " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinX() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinY() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());
            rows.get(firstRow).get(randomContact).getContactRecord().incrementCount(-1);
            rows.get(firstRow).set(randomContact, rows.get(firstRow).get(rows.get(firstRow).size()-1));
            //System.out.println("first row swapped " + currentRow + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinX() + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinY() + " " + rows.get(firstRow).get(randomContact).getLink() +  " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinX() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinY() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());

            long partnerRow;
            int partnerLink = rows.get(firstRow).get(randomContact).getLink();
            if (rows.get(firstRow).get(randomContact).getContactRecord().getBinX() == firstRow) {
                partnerRow = (long) rows.get(firstRow).get(randomContact).getContactRecord().getBinY();
            } else {
                partnerRow = (long) rows.get(firstRow).get(randomContact).getContactRecord().getBinX();
            }
            //System.out.println(firstRow + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinY() + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinX() + " " + partnerRow + " " + partnerLink + " " + symmetricRandomContactPosition);
            rows.get(partnerRow).get(partnerLink).setLink(randomContact);
            //System.out.println("partner link updated " + currentRow + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinX() + " " + rows.get(firstRow).get(randomContact).getContactRecord().getBinY() + " " + rows.get(firstRow).get(randomContact).getLink() +  " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinX() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinY() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());
            rows.get(firstRow).remove(rows.get(firstRow).size()-1);
            //System.out.println("first row removed " + currentRow + " " + rows.get(firstRow).get(Math.min(randomContact,rows.get(firstRow).size()-1)).getContactRecord().getBinX() + " " + rows.get(firstRow).get(Math.min(randomContact,rows.get(firstRow).size()-1)).getContactRecord().getBinY() + " " + rows.get(firstRow).get(Math.min(randomContact,rows.get(firstRow).size()-1)).getLink() +  " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinX() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getContactRecord().getBinY() + " " + rows.get(firstRow).get(rows.get(firstRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());

            rowsums.addTo(firstRow, -1);
            float newFirstRowSum = rowsums.get(firstRow);
            //System.out.println(currentRowSum + " " + rowsums.get(currentRow) + " " + sortedRowSums.get(sortedRowSums.size()-1) );
            rowsumIndices.get(currentRowSum).remove(currentRow);
            Map<Long,Integer> newRow1 = rowsumIndices.get(newFirstRowSum);
            if (newRow1 == null) {
                rowsumIndices.put(newFirstRowSum, new HashMap<>());
                newRow1 = rowsumIndices.get(newFirstRowSum);
                newRow1.put(firstRow,1);
            } else {
                newRow1.put(firstRow,1);
            }
            int switchPlace1 = SumMap.get(firstRowSum);
            float potentialSwitchSum1 = sortedRowSums.get(switchPlace1);
            if (switchPlace1 == sortedRowSums.size()-1) {
                sortedRowSums.set(switchPlace1, newFirstRowSum);
                if (SumMap.get(newFirstRowSum) == null) {
                    SumMap.put(newFirstRowSum, switchPlace1);
                }
            } else {
                sortedRowSums.set(sortedRowSums.size()-1, potentialSwitchSum1);
                sortedRowSums.set(switchPlace1, newFirstRowSum);
                SumMap.put(firstRowSum, switchPlace1+1);
                if (SumMap.get(newFirstRowSum) == null) {
                    SumMap.put(newFirstRowSum, switchPlace1);
                }
            }
            //System.out.println(currentRowSum + " " + rowsums.get(currentRow) + " " + sortedRowSums.get(sortedRowSums.size()-1) );

            if (firstRow != secondRow) {
                float secondRowSum = rowsums.get(secondRow);
                //removedContacts.get(secondRow).add(symmetricRandomContactPosition);
                //remainingContacts.get(secondRow).remove(symmetricRandomContactPosition);
                lastRow = rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinX() == secondRow? rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinY() : rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinX();
                lastLink = rows.get(secondRow).get(rows.get(secondRow).size()-1).getLink();
                //System.out.println("second row last check " + currentRow + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinX() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinY() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());

                rows.get(secondRow).set(symmetricRandomContactPosition, rows.get(secondRow).get(rows.get(secondRow).size()-1));
                //System.out.println("second row swapped " + currentRow + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinX() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinY() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());

                partnerLink = rows.get(secondRow).get(symmetricRandomContactPosition).getLink();
                if (rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() == secondRow) {
                    partnerRow = (long) rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY();
                } else {
                    partnerRow = (long) rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX();
                }
                if (symmetricRandomContactPosition!=rows.get(secondRow).size()-1) {
                    rows.get(partnerRow).get(partnerLink).setLink(symmetricRandomContactPosition);
                }
                //System.out.println("partner link updated " + currentRow + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinX() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getContactRecord().getBinY() + " " + rows.get(secondRow).get(symmetricRandomContactPosition).getLink() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinX() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getContactRecord().getBinY() + " " + rows.get(secondRow).get(rows.get(secondRow).size()-1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());

                rows.get(secondRow).remove(rows.get(secondRow).size()-1);
                if (rows.get(secondRow).size()>0) {
                    //System.out.println("second row removed " + currentRow + " " + rows.get(secondRow).get(Math.min(symmetricRandomContactPosition, rows.get(secondRow).size() - 1)).getContactRecord().getBinX() + " " + rows.get(secondRow).get(Math.min(symmetricRandomContactPosition, rows.get(secondRow).size() - 1)).getContactRecord().getBinY() + " " + rows.get(secondRow).get(Math.min(symmetricRandomContactPosition, rows.get(secondRow).size() - 1)).getLink() + " " + rows.get(secondRow).get(rows.get(secondRow).size() - 1).getContactRecord().getBinX() + " " + rows.get(secondRow).get(rows.get(secondRow).size() - 1).getContactRecord().getBinY() + " " + rows.get(secondRow).get(rows.get(secondRow).size() - 1).getLink() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinX() + " " + rows.get(lastRow).get(lastLink).getContactRecord().getBinY() + " " + rows.get(lastRow).get(lastLink).getLink());
                }
                rowsums.addTo(secondRow, -1);
                //System.out.println(currentRowSum + " " + rowsums.get(currentRow) + " " + sortedRowSums.get(sortedRowSums.size()-1) + " " + secondRowSum );
                float newSecondRowSum = rowsums.get(secondRow);
                //int removeIndex2 = Collections.binarySearch(rowsumIndices.get(secondRowSum),secondRow);
                rowsumIndices.get(secondRowSum).remove(secondRow);
                Map<Long,Integer> newRow2 = rowsumIndices.get(newSecondRowSum);
                if (newRow2 == null) {
                    rowsumIndices.put(newSecondRowSum, new HashMap<>());
                    newRow2 = rowsumIndices.get(newSecondRowSum);
                    newRow2.put(secondRow,1);
                } else {
                    newRow2.put(secondRow,1);
                }
                int switchPlace2 = SumMap.get(secondRowSum);
                if (switchPlace2 == 0) {
                    sortedRowSums.set(0,newSecondRowSum);
                    if (SumMap.get(newSecondRowSum) == null) {
                        SumMap.put(newSecondRowSum,0);
                    }
                    if (sortedRowSums.get(1)==secondRowSum) {
                        SumMap.put(secondRowSum,1);
                    } else {
                        SumMap.remove(secondRowSum);
                    }
                } else {
                    sortedRowSums.set(switchPlace2, newSecondRowSum);
                    if (SumMap.get(newSecondRowSum) == null) {
                        SumMap.put(newSecondRowSum,switchPlace2);
                    }
                    if (sortedRowSums.size() == switchPlace2+1) {
                        SumMap.remove(secondRowSum);
                    } else if (sortedRowSums.get(switchPlace2+1) == secondRowSum) {
                        SumMap.put(secondRowSum, switchPlace2+1);
                    } else {
                        SumMap.remove(secondRowSum);
                    }
                }
            }
            //System.out.println(currentRowSum + " " + rowsums.get(currentRow) + " " + sortedRowSums.get(sortedRowSums.size()-1) );
            currentRows = rowsumIndices.get(sortedRowSums.get(sortedRowSums.size()-1));
            currentRow = currentRows.keySet().iterator().next();
            double maxRatio = (rowsums.get(currentRow)*1.0d)/rowSumThreshold;
            if (thresholds.size()>0 && maxRatio < thresholds.get(thresholds.size()-1)) {
                double passedThreshold = thresholds.get(thresholds.size()-1);
                while (thresholds.size()>0 && maxRatio < thresholds.get(thresholds.size()-1)) {
                    passedThreshold = thresholds.remove(thresholds.size()-1);
                }
                System.out.println("passed threshold: " + passedThreshold + " for matrix: " + chr1 + "-" + chr2 + "(current max sum: " + rowsums.get(currentRow) + " , sum threshold: " + rowSumThreshold + ")");
            }
            //Instant D = Instant.now();
            //System.err.println(Duration.between(C,D).toMillis());
        }

        return contactRecords;

    }

    public int getRandomWithExclusion(Random rnd, int end, List<Integer> exclude) {
        int random = 0;
        try {
            random = rnd.nextInt(end - exclude.size());
        } catch (Exception e) {
            System.err.println(end + " " + exclude.size());
            e.printStackTrace();
        }
        for (int ex : exclude) {
            if (random < ex) {
                break;
            }
            random++;
        }
        return random;
    }*/
}