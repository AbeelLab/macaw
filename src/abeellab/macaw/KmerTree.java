package abeellab.macaw;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;

import java.util.LinkedList;

public class KmerTree {
		
	private KmerTree a = null;
	private KmerTree t = null;
	private KmerTree c = null;
	private KmerTree g = null;


	public KmerTree() {
		super();
	}
	// TODO add exception in case not a/t/c/g/A/T/C/G

	private KmerTree(char base, LinkedList<Character> remaining, int depth, int origin) {
		if (base == 'A' || base == 'a') {
			this.a = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, origin);				
		}
		else if (base == 'T' || base == 't') {
			this.t = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, origin);				
		}
		else if (base == 'C' || base == 'c') {
			this.c = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, origin);				
		}
		else if (base == 'G' || base == 'g') {
			this.g = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, origin);				
		}
	}

	private static KmerTree newNode(char base, LinkedList<Character> remaining, int depth, int origin) {
		//	public static Node newNode(int depth) {
		if (depth > 0) {
			return new KmerTree(base, remaining, depth, origin);
			//			return new Node(depth);
		}
		else {
			//			return new Leaf(base, origin);
			return new Leaf(origin);
		}
	}

	
	public void fill(char[]bases,int position){
		LinkedList<Character> ll=new LinkedList<Character>();
		for(char c:bases)
			ll.add(c);
		fill(ll.removeFirst(),ll,bases.length-1,position);
	}
	
	protected synchronized void fill(char base, LinkedList<Character> remaining, int depth, int position) {
		if (depth > 0) {
			if (base == 'A' || base == 'a') {
				if (this.a == null) {
					this.a = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, position);									
				}
				else {
					this.a.fill(remaining.removeFirst(), remaining, depth-1, position);
				}
			}
			else if (base == 'T' || base == 't') {
				if (this.t == null) {
					this.t = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, position);									
				}
				else {
					this.t.fill(remaining.removeFirst(), remaining, depth-1, position);
				}
			}
			else if (base == 'C' || base == 'c') {
				if (this.c == null) {
					this.c = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, position);
				}
				else {
					this.c.fill(remaining.removeFirst(), remaining, depth-1, position);
				}
			}
			else if (base == 'G' || base == 'g') {
				if (this.g == null) {
					this.g = KmerTree.newNode(remaining.removeFirst(), remaining, depth-1, position);									
				}
				else {
					this.g.fill(remaining.removeFirst(), remaining, depth-1, position);
				}
			}
			else if (base == '\u0000') {
				this.fill(remaining.removeFirst(), remaining, depth, position); // ?????? TODO why???
			}
		}
		else {
			if (this.a != null) {
				this.a.fill(base, remaining, depth, position);
			}
			if (this.t != null) {
				this.t.fill(base, remaining, depth, position);
			}
			if (this.c != null) {
				this.c.fill(base, remaining, depth, position);
			}
			if (this.g != null) {
				this.g.fill(base, remaining, depth, position);
			}
		}
	}


//	public void calculate_distances(Stats s, int origin, int depth) {
//		// TODO Auto-generated method stub
//
//		if (this.a != null) {
//			this.a.calculate_distances(s, origin, depth-1);
//		}
//		if (this.t != null) {
//			this.t.calculate_distances(s, origin, depth-1);
//		}
//		if (this.c != null) {
//			this.c.calculate_distances(s, origin, depth-1);
//		}
//		if (this.g != null) {
//			this.g.calculate_distances(s, origin, depth-1);
//		}
//	}
//
	public String toString() {
		return "A:"+this.a.toString() +"\nT:"+ this.t.toString() +"\nC:"+ this.c.toString() +"\nG:"+ this.g.toString();
	}


}


class Leaf extends KmerTree {

	private TIntArrayList positions = new TIntArrayList(); 


	public Leaf(int position) {
		this.positions.add(position);
	}


	public synchronized void fill(char base, LinkedList<Character> remaining, int depth, int origin) {
			if (!this.positions.contains(origin)) {
				this.positions.add(origin);
			}
	}

}


