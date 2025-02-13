import java.util.ArrayList;
import java.util.HashMap;

//import gurobi.GRBException;
//import gurobi.GRBModel;

//import gurobi.GRBModel;

/*import com.gurobi.gurobi.GRB;
import com.gurobi.gurobi.GRB.DoubleAttr;
import com.gurobi.gurobi.GRBConstr;
import com.gurobi.gurobi.GRBEnv;
import com.gurobi.gurobi.GRBException;
import com.gurobi.gurobi.GRBLinExpr;
import com.gurobi.gurobi.GRBModel;
import com.gurobi.gurobi.GRBVar;*/

public class RMP {
	// the model corresponding to this problem:
	// public GRBModel model;
	// the set of uncovered (also eligible for A, B, and C) vertices:
	public boolean[] uncovered;
	// the covered vertices from A_0:
	public boolean[] covered_A0;
	// the covered vertices from BC_0:
	public boolean[] covered_B0;
	// the covered vertices from BC_0:
	public boolean[] covered_C0;
	// the eligible vertices for A:
	public boolean[] eligible_A;
	// the eligible vertices for B:
	public boolean[] eligible_B;
	// the number of A_0 covered vertices:
	public int nr_A0;
	// the number of BC_0 covered vertices:
	public int nr_B0;
	// the number of A_0 covered vertices:
	public int nr_C0;
	// the number of A eligible vertices:
	public int nr_eligible_A;
	// the number of B eligible vertices:
	public int nr_eligible_B;
	// the index of variables for the active A subsets:
	public HashMap<Integer, VarSet> variables_A;
	// the index of variables for the active B subsets:
	public HashMap<Integer, VarSet> variables_B;
	// the index of variables for the active C subsets:
	public HashMap<Integer, VarSet> variables_C;
	// vertex_A[i] keeps the id(s) for the active A sets containing vertex i:
	public ArrayList<Integer>[] vertex_A;
	// vertex_B[i] keeps the id(s) for the active B sets containing vertex i:
	public ArrayList<Integer>[] vertex_B;
	// vertex_C[i] keeps the id(s) for the active C sets containing vertex i:
	public ArrayList<Integer>[] vertex_C;

	public RMP(boolean[] _uncovered, boolean[] _covered_A0, boolean[] _covered_B0, boolean[] _covered_C0,
			boolean[] _eligible_A, boolean[] _eligible_B, int _nr_A0, int _nr_B0, int _nr_C0, int _nr_eligible_A,
			int _nr_eligible_B, HashMap<Integer, VarSet> _variables_A, HashMap<Integer, VarSet> _variables_B,
			HashMap<Integer, VarSet> _variables_C, ArrayList<Integer>[] _vertex_A, ArrayList<Integer>[] _vertex_B,
			ArrayList<Integer>[] _vertex_C) {

		this.uncovered = _uncovered;
		this.covered_A0 = _covered_A0;
		this.covered_B0 = _covered_B0;
		this.covered_C0 = _covered_C0;
		this.eligible_A = _eligible_A;
		this.eligible_B = _eligible_B;
		this.nr_A0 = _nr_A0;
		this.nr_B0 = _nr_B0;
		this.nr_C0 = _nr_C0;
		this.nr_eligible_A = _nr_eligible_A;
		this.nr_eligible_B = _nr_eligible_B;
		this.variables_A = _variables_A;
		this.variables_B = _variables_B;
		this.variables_C = _variables_C;
		this.vertex_A = _vertex_A;
		this.vertex_B = _vertex_B;
		this.vertex_C = _vertex_C;
	}
}
