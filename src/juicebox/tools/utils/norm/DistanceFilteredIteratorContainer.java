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

import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.iterators.BigContactRecordList;
import javastraw.reader.iterators.IteratorContainer;
import javastraw.reader.iterators.ListIteratorContainer;
import javastraw.reader.iterators.ListOfListIteratorContainer;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class DistanceFilteredIteratorContainer extends IteratorContainer {
    private static boolean USE_DIST_FILTERING = false;
    private static float FILTER_DIST_CUTOFF = 0;
    private final IteratorContainer listContainer;

    public DistanceFilteredIteratorContainer(IteratorContainer ic0, int resolution) {
        super(ic0.getMatrixSize());

        List<List<ContactRecord>> contacts = filterContactRecordsByDistance(ic0, resolution);

        if (contacts.size() == 1) {
            listContainer = new ListIteratorContainer(contacts.get(0), getMatrixSize());
        } else {
            BigContactRecordList bg = BigContactRecordList.populateFromListOfLists(contacts);
            listContainer = new ListOfListIteratorContainer(bg, getMatrixSize());
        }
    }

    public static void setFilterDistance(int distance) {
        FILTER_DIST_CUTOFF = distance;
        USE_DIST_FILTERING = distance > 1;
    }

    public static boolean getUseFilterDistance() {
        return USE_DIST_FILTERING;
    }

    private List<List<ContactRecord>> filterContactRecordsByDistance(IteratorContainer initialIC, int resolution) {
        List<List<ContactRecord>> allRecords = new ArrayList<>();
        List<ContactRecord> tempList = new ArrayList<>();
        int distanceCutoff = (int) Math.ceil(FILTER_DIST_CUTOFF / resolution) + 1;
        int counter = 0;
        Iterator<ContactRecord> iterator = initialIC.getNewContactRecordIterator();
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            if (passesFilter(record, distanceCutoff)) {
                tempList.add(record);
                counter++;
            }
            if (counter > BigContactRecordList.MAX_LIMIT) {
                allRecords.add(tempList);
                tempList = new ArrayList<>();
                counter = 0;
            }
        }
        if (tempList.size() > 0) {
            allRecords.add(tempList);
        }
        return allRecords;
    }

    private boolean passesFilter(ContactRecord record, int distanceCutoff) {
        return Math.abs(record.getBinX() - record.getBinY()) <= distanceCutoff;
    }

    @Override
    public Iterator<ContactRecord> getNewContactRecordIterator() {
        return listContainer.getNewContactRecordIterator();
    }

    @Override
    public ListOfFloatArrays sparseMultiply(ListOfFloatArrays listOfFloatArrays, long l) {
        return listContainer.sparseMultiply(listOfFloatArrays, l);
    }

    @Override
    public void clear() {
        listContainer.clear();
    }
}
