import java.util.ArrayList;
import java.util.Iterator;

public class VarSet {

	public ArrayList<Integer> vertices;
	public double cost;
	public int size;
	// public int id;

	public VarSet() {
		this.vertices = new ArrayList<Integer>();
		this.cost = -1;
		this.size = 0;
		// this.id = -1;
	}

	public boolean check(VarSet set, int n) {
		int i;
		boolean[] currentSet = new boolean[n], newSet = new boolean[n];
		Iterator<Integer> iter1;
		if (this.vertices != null) {
			iter1 = this.vertices.iterator();
			while (iter1.hasNext()) {
				i = iter1.next().intValue();
				currentSet[i] = true;
			}
		}
		if (set.vertices != null) {
			iter1 = set.vertices.iterator();
			while (iter1.hasNext()) {
				i = iter1.next().intValue();
				newSet[i] = true;
			}
		}
		for (i = 0; i < n; i++)
			if (currentSet[i] != newSet[i])
				return false;

		return true;
	}

	public VarSet(VarSet set) {
		this.vertices = new ArrayList<Integer>();
		this.cost = -1;
		this.size = set.vertices.size();
		int i;

		Iterator<Integer> iterator = set.vertices.iterator();
		if (set != null)
			while (iterator.hasNext()) {
				i = iterator.next().intValue();
				this.vertices.add(i);
			}
	}

	public VarSet(boolean[] characteristicSet) {
		this.vertices = new ArrayList<Integer>();
		this.cost = -1;
		int n = characteristicSet.length;
		for (int h = 0; h < n; h++)
			if (characteristicSet[h])
				this.vertices.add(h);
		this.size = this.vertices.size();
	}

	public boolean isEqual(VarSet set, int n) {
		boolean result = true;
		boolean[] charactThis = new boolean[n];
		boolean[] charactThat = new boolean[n];
		int i;
		Iterator<Integer> iter = this.vertices.iterator();
		while (iter.hasNext()) {
			i = iter.next().intValue();
			charactThis[i] = true;
		}
		iter = set.vertices.iterator();
		while (iter.hasNext()) {
			i = iter.next().intValue();
			charactThat[i] = true;
		}
		for (int j = 0; j < n; j++) {
			if (charactThis[j] != charactThat[j])
				return false;
		}
		return result;
	}

	public String toString(int id) {
		String toStr = "";
		Iterator<Integer> iter;
		iter = this.vertices.iterator();
		while (iter.hasNext()) {
			toStr += " " + iter.next().intValue();
		}
		System.out.print("Id: " + id + ", size = " + this.vertices.size() + " | " + toStr);
		System.out.println();
		return toStr;
	}

	public String toString() {
		String toStr = "";
		Iterator<Integer> iter;
		iter = this.vertices.iterator();
		while (iter.hasNext()) {
			toStr += " " + iter.next().intValue();
		}
		System.out.print("size = " + this.vertices.size() + " | " + toStr);
		System.out.println();
		return toStr;
	}

}
