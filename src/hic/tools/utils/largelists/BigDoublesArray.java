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

package hic.tools.utils.largelists;

import com.google.common.util.concurrent.AtomicDouble;
import hic.HiCGlobals;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.tools.ParallelizationTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * can't use <T> because we need to instantiate the array, otherwise that would have been nice
 */
public class BigDoublesArray {

	final long DEFAULT_LENGTH = BigShortsArray.DEFAULT_LENGTH;
	final long overallLength;
	final List<double[]> internalList = new ArrayList<>();

	public BigDoublesArray(long length) {
		this.overallLength = length;
		long tempLength = length;
		while (tempLength > 0) {
			if (tempLength < DEFAULT_LENGTH) {
				internalList.add(new double[(int) tempLength]);
				break;
			} else {
				internalList.add(new double[(int) DEFAULT_LENGTH]);
				tempLength -= DEFAULT_LENGTH;
			}
		}
	}

	public void clear() {
		internalList.clear();
	}

	public double get(long index) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			return internalList.get(pseudoRow)[pseudoCol];
		} else {
			System.err.println("long index exceeds max size of list of arrays while getting: " + index + " " + overallLength);
			Exception ioe = new Exception();
			ioe.printStackTrace();
			return Double.NaN;
		}
	}

	public void set(long index, double value) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			internalList.get(pseudoRow)[pseudoCol] = value;
		} else {
			System.err.println("long index exceeds max size of list of arrays while setting");
		}
	}

	public long getLength() {
		return overallLength;
	}

	public void addTo(long index, double value) {
		if (index < overallLength) {
			int pseudoRow = (int) (index / DEFAULT_LENGTH);
			int pseudoCol = (int) (index % DEFAULT_LENGTH);
			try {
				internalList.get(pseudoRow)[pseudoCol] += value;
			} catch (Exception e) {
				System.err.println(index + " " + pseudoCol);
				e.printStackTrace();
			}
		} else {
			System.err.println("long index exceeds max size of list of arrays while adding: " + index + " " + overallLength);
			Exception ioe = new Exception();
			ioe.printStackTrace();
		}
	}

	public void addValuesFrom(BigDoublesArray other) {
		if (overallLength == other.overallLength) {
			for (int i = 0; i < internalList.size(); i++) {
				for (int j = 0; j < internalList.get(i).length; j++) {
					internalList.get(i)[j] += other.internalList.get(i)[j];
				}
			}
		} else {
			System.err.println("Adding objects of different sizes!");
		}
	}

	public List<double[]> getValues() {
		return internalList;
	}

	public BigFloatsArray convertToFloats() {
		BigFloatsArray newList = new BigFloatsArray(overallLength);
		for (int j = 0; j < internalList.size(); j++) {
			float[] dest = newList.getValues().get(j);
			double[] src = internalList.get(j);
			for (int k = 0; k < dest.length; k++) {
				dest[k] = (float) src[k];
			}
		}
		return newList;
	}

	private int getNumThreads() {
		return Math.min(HiCGlobals.normThreads, internalList.size());
	}

	public static double parCalculateError(BigDoublesArray col, BigDoublesArray scale, BigShortsArray target, BigShortsArray bad) {
		AtomicDouble atomicDouble = new AtomicDouble(0);
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(col.getNumThreads(), () -> {
			int i = index.getAndIncrement();
			double err = 0;
			while (i < col.internalList.size()) {
				double[] colA = col.internalList.get(i);
				double[] scaleA = scale.internalList.get(i);
				short[] targetA = target.internalList.get(i);
				short[] badA = bad.internalList.get(i);

				for (int z = 0; z < colA.length; z++) {
					if (badA[z] == 1) continue;
					double tempErr = Math.abs((colA[z] * scaleA[z] - targetA[z]));
					if (tempErr > err) {
						err = tempErr;
					}
				}
				i = index.getAndIncrement();
			}
			synchronized (atomicDouble) {
				if (err > atomicDouble.get()) {
					atomicDouble.set(err);
				}
			}
		});
		return atomicDouble.get();
	}

	public static double calculateError90(BigDoublesArray col, BigDoublesArray scale,
										  BigShortsArray target, BigShortsArray bad) {
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for (int i = 0; i < col.internalList.size(); i++) {
			double[] colA = col.internalList.get(i);
			double[] scaleA = scale.internalList.get(i);
			short[] targetA = target.internalList.get(i);
			short[] badA = bad.internalList.get(i);

			for (int z = 0; z < colA.length; z++) {
				if (badA[z] == 1) continue;
				double tempErr = Math.abs((colA[z] * scaleA[z] - targetA[z]));
				stats.addValue(tempErr);
			}
		}
		return stats.getPercentile(90);
	}

	public static double parCalculateConvergenceError(BigDoublesArray calculatedVectorB, BigDoublesArray current,
													  BigShortsArray bad) {
		AtomicDouble atomicDouble = new AtomicDouble(0);
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(current.getNumThreads(), () -> {
			int i = index.getAndIncrement();
			double err = 0;
			while (i < current.internalList.size()) {

				double[] calcA = calculatedVectorB.internalList.get(i);
				double[] currA = current.internalList.get(i);
				short[] badA = bad.internalList.get(i);

				for (int z = 0; z < calcA.length; z++) {
					if (badA[z] == 1) continue;
					double tempErr = Math.abs(calcA[z] - currA[z]);
					if (tempErr > err) {
						err = tempErr;
					}
				}
				i = index.getAndIncrement();
			}
			synchronized (atomicDouble) {
				if (err > atomicDouble.get()) {
					atomicDouble.set(err);
				}
			}
		});
		return atomicDouble.get();
	}

	public BigDoublesArray deepClone() {
		BigDoublesArray clone = new BigDoublesArray(overallLength);
		for (int k = 0; k < internalList.size(); k++) {
			System.arraycopy(internalList.get(k), 0, clone.internalList.get(k), 0, internalList.get(k).length);
		}
		return clone;
	}

	public BigShortsArray deepCovertedClone() {
		BigShortsArray clone = new BigShortsArray(overallLength);
		for (long p = 0; p < overallLength; p++) {
			clone.set(p, (short) get(p));
		}
		return clone;
	}

	public void parSetToGeoMean(BigDoublesArray a, BigDoublesArray b) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				double[] result = internalList.get(i);
				double[] a1 = a.internalList.get(i);
				double[] b1 = b.internalList.get(i);
				for (int p = 0; p < result.length; p++) {
					result[p] = (float) Math.sqrt(a1[p] * b1[p]);
				}
				i = index.getAndIncrement();
			}
		});
	}

	public void parSetTo(BigDoublesArray srcArrays) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				double[] dest = internalList.get(i);
				double[] src = srcArrays.internalList.get(i);
				System.arraycopy(src, 0, dest, 0, dest.length);
				i = index.getAndIncrement();
			}
		});
	}

	public void parMultiplyByOneMinus(BigShortsArray array) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				double[] orig = internalList.get(i);
				short[] arr = array.internalList.get(i);
				for (int p = 0; p < orig.length; p++) {
					orig[p] *= (1 - arr[p]);
				}
				i = index.getAndIncrement();
			}
		});
	}

	public void parMultiplyBy(BigDoublesArray dv) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				double[] orig = internalList.get(i);
				double[] arr = dv.internalList.get(i);
				for (int p = 0; p < orig.length; p++) {
					orig[p] *= arr[p];
				}
				i = index.getAndIncrement();
			}
		});
	}

	public void parSetToDivision(BigShortsArray num, BigDoublesArray denom) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				double[] orig = internalList.get(i);
				short[] num1 = num.internalList.get(i);
				double[] denom1 = denom.internalList.get(i);
				for (int p = 0; p < orig.length; p++) {
					orig[p] = num1[p] / denom1[p];
				}
				i = index.getAndIncrement();
			}
		});
	}

	public void parScaleByRatio(BigShortsArray num, BigDoublesArray denom) {
		AtomicInteger index = new AtomicInteger();
		ParallelizationTools.launchParallelizedCode(getNumThreads(), () -> {
			int i = index.getAndIncrement();
			while (i < internalList.size()) {
				double[] orig = internalList.get(i);
				short[] num1 = num.internalList.get(i);
				double[] denom1 = denom.internalList.get(i);
				for (int p = 0; p < orig.length; p++) {
					orig[p] *= (num1[p] / denom1[p]);
				}
				i = index.getAndIncrement();
			}
		});
	}

	public ListOfFloatArrays convertToListOfFloat() {
		ListOfFloatArrays clone = new ListOfFloatArrays(overallLength);
		for (long k = 0; k < getLength(); k++) {
			clone.set(k, (float) get(k));
		}
		return clone;
	}
}
