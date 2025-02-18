import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Stack;

import org.graph4j.Digraph;
import org.graph4j.GraphBuilder;
import org.graph4j.alg.flow.PushRelabelMaximumFlow;

//import gurobi.GRB;
//import gurobi.GRB.DoubleAttr;
//import gurobi.GRBConstr;
//import gurobi.GRBEnv;
//import gurobi.GRBException;
//import gurobi.GRBLinExpr;
//import gurobi.GRBModel;
//import gurobi.GRBVar;

import com.gurobi.gurobi.GRB;
import com.gurobi.gurobi.GRB.DoubleAttr;
import com.gurobi.gurobi.GRBConstr;
import com.gurobi.gurobi.GRBEnv;
import com.gurobi.gurobi.GRBException;
import com.gurobi.gurobi.GRBLinExpr;
import com.gurobi.gurobi.GRBModel;
import com.gurobi.gurobi.GRBVar;

import java.util.Map.Entry;
import java.util.concurrent.TimeUnit;

//TODO 1: <n &&..
//TODO 2: don't add C-subsets if |C|>= best_nr_C

public class BranchPriceImp {
	// the order and the size of G, the upper bound for A and B, the max id for a
	// variable, best card(C):
	public static int n, m, b, varId = 0, bestSeparatorSize, best_nr_C = n, totalNoOfVars = 0, rootNoOfVars,
			init_nr_C = 0, time_max = 3600;
	// best solution: 1 for A, 2 for B, 3 for C:
	public static int[] bestSolution;

	// the adjacency matrix:
	public static boolean[][] adjacency;
	public static double alpha_star, rootTime, rootObj;
	// the array of degrees:
	public static int[] degrees;
	// the stack of RMPs:
	public static Stack<RMP> stackRMP = new Stack<RMP>();
	// the reduced cost precision:
	public static double precisionRC = 1e-8;
	// the subproblem weights precision:
	public static double precisionWeight = 1e-9;
	// variables precision:
	public static double precisionVar = 1e-9;
	// optimal objective precision:
	public static double precisionOpt = 1e-8;
	// VarSet weight precision:
	public static double precisionSetWeight = 1e-8;

	public static boolean correctSolution = false;

	public static String fileName = "", fileNameInstances = "";// to be initialized
	public static String path = "../data/VSP/";
	public static String path_results = "../data/VSP/results/";
	public static String path_instances = "../data/VSP/instances/";
	public static String path_solutions = "../data/VSP/solutions/";

	public static void createDir(String path) {
		File file = new File(path);
		boolean bool = file.mkdir();
		if (bool) {
			System.out.println("Results directory created successfully");
		} else {
			System.out.println("Sorry couldn't create results directory. Maybe already there?");
		}
	}

	public static int[] readGraphSize(String dataFile) {
		int[] size = new int[2];
		String text = null;
		File file = new File(path_instances + dataFile);
		BufferedReader reader = null;
		String[] nors = null;

		try {
			reader = new BufferedReader(new FileReader(file));
			while ((text = reader.readLine()) != null && text.startsWith("c")) {
			}
			if (text.startsWith("p")) {
				nors = text.split("\\s+", 0);
				size[0] = Integer.valueOf(nors[2]);
				size[1] = Integer.valueOf(nors[3]);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return size;
	}

	public static int[] readGraphSizeTable(String dataFile) {
		int[] size = new int[2];
		String text = null;
		File file = new File(path_instances + dataFile);
		BufferedReader reader = null;
		String[] nors = null;

		try {
			reader = new BufferedReader(new FileReader(file));
			text = reader.readLine();
			nors = text.split("\\s+", 0);
			size[0] = Integer.parseInt(nors[1]);
			size[1] = Integer.parseInt(nors[2]);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return size;
	}

	public static int readGraphFromFile(String dataFile) {
		int[] size = readGraphSize(dataFile);
		m = size[1];
		n = size[0];
		int max = 0;
		degrees = new int[n];
		adjacency = new boolean[n][n];
		int i, j;
		String text = null;
		File file = new File(path_instances + dataFile);
		BufferedReader reader = null;
		String[] nors = null;
		System.out.println("Input file: " + dataFile + ", " + n + " vertices" + ", " + m + " edges ("
				+ (n * (n - 1) / 2 - m) + " non-adjacencies)");
		try {
			reader = new BufferedReader(new FileReader(file));
			while ((text = reader.readLine()) != null && text.startsWith("c")) {
			}
			for (int k = 0; k < m; k++) {
				if ((text = reader.readLine()) != null && text.startsWith("e")) {
					nors = text.split("\\s+", 0);
					i = Integer.valueOf(nors[1]);
					j = Integer.valueOf(nors[2]);
					adjacency[i - 1][j - 1] = true;
					adjacency[j - 1][i - 1] = true;
					degrees[i - 1]++;
				}
			}
			for (int k = 0; k < n; k++)
				if (max < degrees[k])
					max = degrees[k];

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return max;
	}

	public static int readGraphFromFileTable(String dataFile) {
		int[] size = readGraphSizeTable(dataFile);
		m = size[1];
		n = size[0];
		int max = 0;
		degrees = new int[n];
		adjacency = new boolean[n][n];
		int i, j;
		String text = null;
		File file = new File(path_instances + dataFile);
		BufferedReader reader = null;
		String[] nors = null;
		System.out.println("Input file: " + dataFile + ", " + n + " vertices" + ", " + m + " edges ("
				+ (n * (n - 1) / 2 - m) + " non-adjacencies)");
		try {
			reader = new BufferedReader(new FileReader(file));
			while ((text = reader.readLine()) != null && !text.startsWith("e")) {
			}
			for (int k = 0; k < m; k++) {
				if ((text = reader.readLine()) != null && text.startsWith("e")) {
					nors = text.split("\\s+", 0);
					i = Integer.parseInt(nors[1]);
					j = Integer.parseInt(nors[2]);
					adjacency[i][j] = true;
					adjacency[j][i] = true;
					degrees[i]++;
				}
			}
			for (int k = 0; k < n; k++)
				if (adjacency[n - 1][k])
					degrees[n - 1]++;

			for (int k = 0; k < n; k++)
				if (max < degrees[k])
					max = degrees[k];

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return max;
	}

	public static void createFile(String fileName) {
		try {
			File file = new File(fileName);
			if (file.createNewFile()) {
				System.out.println("File created: " + file.getName());
			} else {
				System.out.println("File already exists.");
			}
		} catch (IOException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}
	}

	public static ArrayList<String> createListOfInstances(String fileName) {
		ArrayList<String> instances = new ArrayList<String>();
		File file = new File(path + fileName);
		BufferedReader reader = null;
		String[] nors = null;
		String text = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			while ((text = reader.readLine()) != null)
				instances.add(text);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return instances;
	}

	public static void readSolution(RMP rmp, GRBModel model) throws GRBException {
		GRBVar var;
		Iterator<Entry<Integer, VarSet>> iterVarSet;
		double varValue;
		String varName;
		if (rmp.variables_A != null) {
			System.out.println("x variables");
			iterVarSet = rmp.variables_A.entrySet().iterator();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				varName = "x_" + pairE.getKey();
				var = model.getVarByName(varName);
				varValue = var.get(GRB.DoubleAttr.X);// model.getVarByName("x_" + pairE.getKey()).get(GRB.DoubleAttr.X);
				if (varValue > precisionVar)
					System.out.println(varName + " = " + varValue + " ");
			}
		}
		if (rmp.variables_B != null) {
			// System.out.println("y variables");
			iterVarSet = rmp.variables_B.entrySet().iterator();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				varName = "y_" + pairE.getKey();
				var = model.getVarByName(varName);
				varValue = var.get(GRB.DoubleAttr.X);
				if (varValue > precisionVar)
					System.out.print(varName + " = " + varValue + "");
			}
		}
		if (rmp.variables_C != null) {
			System.out.println("z variables");
			iterVarSet = rmp.variables_C.entrySet().iterator();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				varName = "z_" + pairE.getKey();
				var = model.getVarByName(varName);
				varValue = var.get(GRB.DoubleAttr.X);
				if (varValue > precisionVar)
					System.out.print(varName + " = " + varValue + " ");
			}
		}
	}

	public static Digraph buildDigraph() {
		Digraph g = GraphBuilder.numVertices(2 * n).estimatedNumEdges(n + 2 * m).buildDigraph();
		// vertex i is transformed in (2 * i) = a_i and (2 * i + 1) = b_i
		for (int i = 0; i < n; i++) {
			g.addEdge(2 * i, 2 * i + 1, 1);
			for (int j = 0; j < i; j++)
				if (adjacency[i][j]) {
					g.addEdge(2 * i + 1, 2 * j, 2 * n);
					g.addEdge(2 * j + 1, 2 * i, 2 * n);
				}
		}
		return g;
	}

	public static GRBModel buildLPModel(RMP rmp) throws GRBException {
		// build (the initial) LP model for the Vertex Separator Problem (VSP)
		GRBVar var;
		GRBEnv env = new GRBEnv();
		GRBModel model = new GRBModel(env);
		GRBLinExpr expr1, expr2_A = new GRBLinExpr(), expr2_B = new GRBLinExpr(), expr2_C = new GRBLinExpr(),
				expr3 = new GRBLinExpr(), objFunction = new GRBLinExpr();
		Iterator<Entry<Integer, VarSet>> iterVarSet;
		Iterator<Integer> iterVarSetId;
		int id;
		boolean found = false, non_empty_AB = false;

		if (rmp.variables_A != null) {
			non_empty_AB = true;
			iterVarSet = rmp.variables_A.entrySet().iterator();
			expr1 = new GRBLinExpr();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				var = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "x_" + pairE.getKey());
				expr1.addTerm(pairE.getValue().size, var);
				expr2_A.addTerm(1.0, var);
				expr3.addTerm(pairE.getValue().size, var);
			}
			model.addConstr(expr1, GRB.LESS_EQUAL, b - rmp.nr_A0, "constraint_A");
			if (rmp.nr_A0 == 0)
				model.addConstr(expr2_A, GRB.EQUAL, 1, "one_subset_A");
			else
				model.addConstr(expr2_A, GRB.LESS_EQUAL, 1, "one_subset_A");
		}

		if (rmp.variables_B != null) {
			non_empty_AB = true;
			iterVarSet = rmp.variables_B.entrySet().iterator();
			expr1 = new GRBLinExpr();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				var = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "y_" + pairE.getKey());
				expr1.addTerm(pairE.getValue().size, var);
				expr2_B.addTerm(1.0, var);
				expr3.addTerm(-pairE.getValue().size, var);
			}
			model.addConstr(expr1, GRB.LESS_EQUAL, b - rmp.nr_B0, "constraint_B");
			if (rmp.nr_B0 == 0)
				model.addConstr(expr2_B, GRB.EQUAL, 1, "one_subset_B");
			else
				model.addConstr(expr2_B, GRB.LESS_EQUAL, 1, "one_subset_B");
		}

		if (non_empty_AB)
			model.addConstr(expr3, GRB.GREATER_EQUAL, rmp.nr_B0 - rmp.nr_A0, "asymmetric");

		if (rmp.variables_C != null) {
			iterVarSet = rmp.variables_C.entrySet().iterator();
			expr1 = new GRBLinExpr();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				var = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "z_" + pairE.getKey());
				objFunction.addTerm(pairE.getValue().vertices.size(), var);
				expr2_C.addTerm(1.0, var);
			}
			if (rmp.nr_C0 == 0)
				model.addConstr(expr2_C, GRB.EQUAL, 1, "one_subset_C");
			else
				model.addConstr(expr2_C, GRB.LESS_EQUAL, 1, "one_subset_C");
		}

		model.setObjective(objFunction);
		model.update();

		for (int i = 0; i < n; i++)
			if (rmp.uncovered[i]) {
				expr1 = new GRBLinExpr();
				if (rmp.vertex_A != null && rmp.vertex_A[i] != null && rmp.vertex_A[i].size() > 0) {
					iterVarSetId = rmp.vertex_A[i].iterator();
					found = true;
					while (iterVarSetId.hasNext()) {
						id = iterVarSetId.next();
						expr1.addTerm(1.0, model.getVarByName("x_" + id));
					}
				}
				if (rmp.vertex_B != null && rmp.vertex_B[i] != null && rmp.vertex_B[i].size() > 0) {
					iterVarSetId = rmp.vertex_B[i].iterator();
					found = true;
					while (iterVarSetId.hasNext()) {
						id = iterVarSetId.next();
						expr1.addTerm(1.0, model.getVarByName("y_" + id));
					}
				}
				if (rmp.vertex_C != null && rmp.vertex_C[i] != null && rmp.vertex_C[i].size() > 0) {
					iterVarSetId = rmp.vertex_C[i].iterator();
					found = true;
					while (iterVarSetId.hasNext()) {
						id = iterVarSetId.next();
						expr1.addTerm(1.0, model.getVarByName("z_" + id));
					}
				}
				if (found)
					model.addConstr(expr1, GRB.EQUAL, 1.0, "vertex_" + i);
				else
					System.out.println("No constraint for the vertex " + i);
			}

		for (int u = 0; u < n; u++)
			if (rmp.eligible_A[u])
				for (int v = u + 1; v < n; v++)
					if (rmp.eligible_B[v] && adjacency[u][v]) {
						found = false;
						expr1 = new GRBLinExpr();
						if (rmp.vertex_A != null && rmp.vertex_A[u] != null && rmp.vertex_A[u].size() > 0) {
							iterVarSetId = rmp.vertex_A[u].iterator();
							found = true;
							while (iterVarSetId.hasNext()) {
								id = iterVarSetId.next();
								expr1.addTerm(1.0, model.getVarByName("x_" + id));
							}
						}
						if (rmp.vertex_B != null && rmp.vertex_B[v] != null && rmp.vertex_B[v].size() > 0) {
							iterVarSetId = rmp.vertex_B[v].iterator();
							found = true;
							while (iterVarSetId.hasNext()) {
								id = iterVarSetId.next();
								expr1.addTerm(1.0, model.getVarByName("y_" + id));
							}
						}
						if (found)
							model.addConstr(expr1, GRB.LESS_EQUAL, 1.0, "A_B_edge_" + u + "_" + v);// A_B_edge_
					}
		for (int u = 0; u < n; u++)
			if (rmp.eligible_B[u])
				for (int v = u + 1; v < n; v++)
					if (rmp.eligible_A[v] && adjacency[u][v]) {
						found = false;
						expr1 = new GRBLinExpr();
						if (rmp.vertex_A != null && rmp.vertex_A[v] != null && rmp.vertex_A[v].size() > 0) {
							iterVarSetId = rmp.vertex_A[v].iterator();
							found = true;
							while (iterVarSetId.hasNext()) {
								id = iterVarSetId.next();
								expr1.addTerm(1.0, model.getVarByName("x_" + id));
							}
						}
						if (rmp.vertex_B != null && rmp.vertex_B[u] != null && rmp.vertex_B[u].size() > 0) {
							iterVarSetId = rmp.vertex_B[u].iterator();
							found = true;
							while (iterVarSetId.hasNext()) {
								id = iterVarSetId.next();
								expr1.addTerm(1.0, model.getVarByName("y_" + id));
							}
						}
						if (found)
							model.addConstr(expr1, GRB.LESS_EQUAL, 1.0, "B_A_edge_" + u + "_" + v);// B_A_edge_
					}
		model.update();
		model.write(path_results + fileName + "/model_pb.lp");
		env.set(GRB.IntParam.OutputFlag, 0);
		model.set(GRB.IntParam.OutputFlag, 0);
		model.setObjective(objFunction);
		model.set(GRB.IntAttr.ModelSense, GRB.MINIMIZE);
		model.update();

		return model;
	}

	public static int branchOnVertex(GRBModel model, RMP rmp) throws GRBException {
		// "model" was already solved;
		// We test if the current solution is integer;
		// if is not integer, then we choose a vertex u for branching.

		boolean found_A, found_B, found_C;
		int u = -1, id;
		double sum_A = 0, sum_B = 0, sum_C = 0, value;
		Iterator<Integer> iterVarSetId;

		for (int v = 0; v < n; v++)
			if (rmp.uncovered[v]) {
				sum_A = -1;
				found_A = false;
				if (rmp.vertex_A != null && rmp.vertex_A[v] != null) {
					sum_A = 0;
					iterVarSetId = rmp.vertex_A[v].iterator();
					while (iterVarSetId.hasNext()) {
						id = iterVarSetId.next();
						value = model.getVarByName("x_" + id).get(GRB.DoubleAttr.X);
						sum_A += value;
						if (value > precisionVar && value < 1 - precisionVar)
							u = v;
					}
				}
				if (sum_A != -1 && sum_A > precisionVar && sum_A < 1 - precisionVar)
					found_A = true;

				sum_B = -1;
				found_B = false;
				if (rmp.vertex_B != null && rmp.vertex_B[v] != null) {
					sum_B = 0;
					iterVarSetId = rmp.vertex_B[v].iterator();
					while (iterVarSetId.hasNext()) {
						id = iterVarSetId.next();
						value = model.getVarByName("y_" + id).get(GRB.DoubleAttr.X);
						sum_B += value;
						if (value > precisionVar && value < 1 - precisionVar)
							u = v;
					}
				}
				if (sum_B != -1 && sum_B > precisionVar && sum_B < 1 - precisionVar)
					found_B = true;

				sum_C = -1;
				found_C = false;
				if (rmp.vertex_C != null && rmp.vertex_C[v] != null) {
					sum_C = 0;
					iterVarSetId = rmp.vertex_C[v].iterator();
					while (iterVarSetId.hasNext()) {
						id = iterVarSetId.next();
						value = model.getVarByName("z_" + id).get(GRB.DoubleAttr.X);
						sum_C += value;
						if (value > precisionVar && value < 1 - precisionVar)
							u = v;
					}
				}
				if (sum_C != -1 && sum_C > precisionVar && sum_C < 1 - precisionVar)
					found_C = true;

				if (found_A || found_B || found_C)
					return u;
			}
		return u;
	}

	public static void branchNImproved(RMP rmp, int u) {
		// When this method is called u != -1.
		RMP rmp_a, rmp_b, rmp_c;
		VarSet vSet, vSet_A, vSet_B, vSet_C;
		Iterator<Entry<Integer, VarSet>> iterVarSet;
		Iterator<Integer> iter;
		int vId, v, c, lowerBound_AB = n - b - best_nr_C + 1, _nr_A0_a, _nr_A0_b, _nr_A0_c, _nr_B0_a, _nr_B0_b,
				_nr_B0_c, _nr_C0_a, _nr_C0_b, _nr_C0_c, _nr_eligible_A_a, _nr_eligible_A_b, _nr_eligible_A_c,
				_nr_eligible_B_a, _nr_eligible_B_b, _nr_eligible_B_c;
		boolean eligible_a = false, eligible_b = false;
		boolean[] u_neighborhood = new boolean[n];
		int[] newIntSolution = new int[n];

		for (int w = 0; w < n; w++)
			if (adjacency[u][w] && rmp.uncovered[w])
				u_neighborhood[w] = true;
		if (rmp.nr_A0 + rmp.nr_B0 + rmp.nr_C0 <= n - 1) {
			if (rmp.nr_A0 <= b - 1 || rmp.nr_B0 <= b - 1) {
				// C branch:
				boolean[] _uncovered_c = new boolean[n], _covered_A0_c = new boolean[n], _covered_B0_c = new boolean[n],
						_covered_C0_c = new boolean[n], _eligible_A_c = new boolean[n], _eligible_B_c = new boolean[n];

				for (int w = 0; w < n; w++) {
					_uncovered_c[w] = rmp.uncovered[w];
					_covered_A0_c[w] = rmp.covered_A0[w];
					_covered_B0_c[w] = rmp.covered_B0[w];
					_covered_C0_c[w] = rmp.covered_C0[w];
					_eligible_A_c[w] = rmp.eligible_A[w];
					_eligible_B_c[w] = rmp.eligible_B[w];
				}
				_uncovered_c[u] = false;
				_covered_C0_c[u] = true;
				_eligible_A_c[u] = false;
				_eligible_B_c[u] = false;

				_nr_A0_c = rmp.nr_A0;
				_nr_B0_c = rmp.nr_B0;
				_nr_C0_c = rmp.nr_C0 + 1;

				_nr_eligible_A_c = rmp.nr_eligible_A;
				if (rmp.eligible_A[u])
					_nr_eligible_A_c = rmp.nr_eligible_A - 1;

				_nr_eligible_B_c = rmp.nr_eligible_B;
				if (rmp.eligible_B[u])
					_nr_eligible_B_c = rmp.nr_eligible_B - 1;

				for (int w = 0; w < n; w++)
					if (!_eligible_A_c[w] && !_eligible_B_c[w] && _uncovered_c[w]) {
						_uncovered_c[w] = false;
						_covered_C0_c[w] = true;
						_nr_C0_c++;
					}
				if (_nr_C0_c < best_nr_C) {
					if (_nr_eligible_A_c > 0 && _nr_eligible_B_c > 0) {
						HashMap<Integer, VarSet> _variables_A_c = new HashMap<Integer, VarSet>(),
								_variables_B_c = new HashMap<Integer, VarSet>(),
								_variables_C_c = new HashMap<Integer, VarSet>();
						ArrayList<Integer>[] _vertex_A_c = new ArrayList[n], _vertex_B_c = new ArrayList[n],
								_vertex_C_c = new ArrayList[n];

						if (rmp.variables_A != null) {
							for (int w = 0; w < n; w++)
								_vertex_A_c[w] = new ArrayList<Integer>();
							iterVarSet = rmp.variables_A.entrySet().iterator();
							while (iterVarSet.hasNext()) {
								HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
										.next();
								vSet = (VarSet) pairE.getValue();
								vId = pairE.getKey();
								if (vSet != null && vSet.vertices.size() > 0) {
									vSet_A = new VarSet();
									iter = vSet.vertices.iterator();
									while (iter.hasNext()) {
										v = iter.next();
										if (_eligible_A_c[v])
											vSet_A.vertices.add(v);
									}
									vSet_A.size = vSet_A.vertices.size();
									if (vSet_A != null && vSet_A.vertices.size() > 0) {
										_variables_A_c.put(vId, vSet_A);
										iter = vSet_A.vertices.iterator();
										while (iter.hasNext())
											_vertex_A_c[iter.next()].add(vId);
									}
								}
							}
						}
						if (rmp.variables_B != null) {
							for (int w = 0; w < n; w++)
								_vertex_B_c[w] = new ArrayList<Integer>();
							iterVarSet = rmp.variables_B.entrySet().iterator();
							while (iterVarSet.hasNext()) {
								HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
										.next();
								vSet = (VarSet) pairE.getValue();
								vId = pairE.getKey();
								if (vSet != null && vSet.vertices.size() > 0) {
									vSet_B = new VarSet();
									iter = vSet.vertices.iterator();
									while (iter.hasNext()) {
										v = iter.next();
										if (_eligible_B_c[v])
											vSet_B.vertices.add(v);
									}
									vSet_B.size = vSet_B.vertices.size();
									if (vSet_B != null && vSet_B.vertices.size() > 0) {
										_variables_B_c.put(vId, vSet_B);
										iter = vSet_B.vertices.iterator();
										while (iter.hasNext())
											_vertex_B_c[iter.next()].add(vId);
									}
								}
							}
						}
						if (rmp.variables_C != null) {
							for (int w = 0; w < n; w++)
								_vertex_C_c[w] = new ArrayList<Integer>();
							iterVarSet = rmp.variables_C.entrySet().iterator();
							while (iterVarSet.hasNext()) {
								HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
										.next();
								vSet = (VarSet) pairE.getValue();
								vId = pairE.getKey();
								if (vSet != null && vSet.vertices.size() > 0) {
									vSet_C = new VarSet();
									iter = vSet.vertices.iterator();
									while (iter.hasNext()) {
										v = iter.next();
										if (_uncovered_c[v])
											vSet_C.vertices.add(v);
									}
									vSet_C.size = vSet_C.vertices.size();
									if (vSet_C != null && vSet_C.vertices.size() > 0
											&& vSet_C.vertices.size() <= best_nr_C) {
										_variables_C_c.put(vId, vSet_C);
										iter = vSet_C.vertices.iterator();
										while (iter.hasNext())
											_vertex_C_c[iter.next()].add(vId);
									}
								}
							}
						}
						rmp_c = new RMP(_uncovered_c, _covered_A0_c, _covered_B0_c, _covered_C0_c, _eligible_A_c,
								_eligible_B_c, _nr_A0_c, _nr_B0_c, _nr_C0_c, _nr_eligible_A_c, _nr_eligible_B_c,
								_variables_A_c, _variables_B_c, _variables_C_c, _vertex_A_c, _vertex_B_c, _vertex_C_c);
						stackRMP.push(rmp_c);

					} else {// _nr_eligible_A_c == 0 or _nr_eligible_B_c == 0
						if (_nr_eligible_A_c == 0 && _nr_eligible_B_c > 0) {// all the vertices in V' are B eligible;
							if (_nr_C0_c < best_nr_C) {
								newIntSolution = new int[n];
								for (int w = 0; w < n; w++) {
									if (_covered_A0_c[w])
										newIntSolution[w] = 1;
									if (_covered_B0_c[w] || _eligible_B_c[w])
										newIntSolution[w] = 2;
									if (_covered_C0_c[w])
										newIntSolution[w] = 3;
								}
								if (checkIntegerSolution(newIntSolution))
									bestSolution = newIntSolution;
							}
						}
						if (_nr_eligible_A_c > 0 && _nr_eligible_B_c == 0) {// all the vertices in V' are A eligible;
							if (_nr_C0_c < best_nr_C) {
								newIntSolution = new int[n];
								for (int w = 0; w < n; w++) {
									if (_covered_A0_c[w] || _eligible_A_c[w])
										newIntSolution[w] = 1;
									if (_covered_B0_c[w])
										newIntSolution[w] = 2;
									if (_covered_C0_c[w])
										newIntSolution[w] = 3;
								}
								if (checkIntegerSolution(newIntSolution))
									bestSolution = newIntSolution;
							}
						}
						if (_nr_eligible_A_c == 0 && _nr_eligible_B_c == 0) {// there are no remaining eligible
																				// vertices;
							if (_nr_C0_c < best_nr_C) {
								newIntSolution = new int[n];
								for (int w = 0; w < n; w++) {
									if (_covered_A0_c[w])
										newIntSolution[w] = 1;
									if (_covered_B0_c[w])
										newIntSolution[w] = 2;
									if (_covered_C0_c[w])
										newIntSolution[w] = 3;
								}
								if (checkIntegerSolution(newIntSolution))
									bestSolution = newIntSolution;
							}
						}
					}
				}
				// A branch:
				if (rmp.eligible_A[u] && rmp.nr_A0 <= b - 1) {
					eligible_a = true;
					boolean[] _uncovered_a = new boolean[n], _covered_A0_a = new boolean[n],
							_covered_B0_a = new boolean[n], _covered_C0_a = new boolean[n],
							_eligible_A_a = new boolean[n], _eligible_B_a = new boolean[n];

					for (int w = 0; w < n; w++) {
						_uncovered_a[w] = rmp.uncovered[w];
						_covered_A0_a[w] = rmp.covered_A0[w];
						_covered_B0_a[w] = rmp.covered_B0[w];
						_covered_C0_a[w] = rmp.covered_C0[w];
						_eligible_A_a[w] = rmp.eligible_A[w];
						_eligible_B_a[w] = rmp.eligible_B[w];
					}
					_uncovered_a[u] = false;
					_covered_A0_a[u] = true;
					_eligible_A_a[u] = false;

					_nr_A0_a = rmp.nr_A0 + 1;
					_nr_B0_a = rmp.nr_B0;
					_nr_C0_a = rmp.nr_C0;

					_nr_eligible_A_a = rmp.nr_eligible_A - 1;
					_nr_eligible_B_a = rmp.nr_eligible_B;
					if (_eligible_B_a[u])
						_nr_eligible_B_a--;
					_eligible_B_a[u] = false;

					for (int w = 0; w < n; w++)
						if (u_neighborhood[w] && _eligible_B_a[w]) {
							_eligible_B_a[w] = false;
							_nr_eligible_B_a--;
						}
					for (int w = 0; w < n; w++)
						if (!_eligible_A_a[w] && !_eligible_B_a[w] && _uncovered_a[w]) {
							_uncovered_a[w] = false;
							_covered_C0_a[w] = true;
							_nr_C0_a++;
						}
					if (_nr_C0_a < best_nr_C) {
						if (_nr_eligible_A_a > 0 && _nr_eligible_B_a > 0) {
							HashMap<Integer, VarSet> _variables_A_a = new HashMap<Integer, VarSet>(),
									_variables_B_a = new HashMap<Integer, VarSet>(),
									_variables_C_a = new HashMap<Integer, VarSet>();

							ArrayList<Integer>[] _vertex_A_a = new ArrayList[n], _vertex_B_a = new ArrayList[n],
									_vertex_C_a = new ArrayList[n];

							if (rmp.variables_A != null && _nr_eligible_A_a > 0) {
								for (int w = 0; w < n; w++)
									_vertex_A_a[w] = new ArrayList<Integer>();
								iterVarSet = rmp.variables_A.entrySet().iterator();
								while (iterVarSet.hasNext()) {
									HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
											.next();
									vSet = (VarSet) pairE.getValue();
									vId = pairE.getKey();
									if (vSet != null && vSet.vertices.size() >= lowerBound_AB) {
										vSet_A = new VarSet();
										iter = vSet.vertices.iterator();
										while (iter.hasNext()) {
											v = iter.next();
											if (_eligible_A_a[v])
												vSet_A.vertices.add(v);
										}
										vSet_A.size = vSet_A.vertices.size();
										if (vSet_A != null && vSet_A.vertices.size() >= lowerBound_AB
												&& vSet_A.vertices.size() > 0) {
											_variables_A_a.put(vId, vSet_A);
											iter = vSet_A.vertices.iterator();
											while (iter.hasNext())
												_vertex_A_a[iter.next()].add(vId);
										}
									}
								}
							}
							if (rmp.variables_B != null && _nr_eligible_B_a > 0) {
								for (int w = 0; w < n; w++)
									_vertex_B_a[w] = new ArrayList<Integer>();
								iterVarSet = rmp.variables_B.entrySet().iterator();
								while (iterVarSet.hasNext()) {
									HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
											.next();
									vSet = (VarSet) pairE.getValue();
									vId = pairE.getKey();
									if (vSet != null && vSet.vertices.size() >= lowerBound_AB) {
										vSet_B = new VarSet();
										iter = vSet.vertices.iterator();
										while (iter.hasNext()) {
											v = iter.next();
											if (_eligible_B_a[v])
												vSet_B.vertices.add(v);
										}
										vSet_B.size = vSet_B.vertices.size();
										if (vSet_B != null && vSet_B.vertices.size() >= lowerBound_AB
												&& vSet_B.vertices.size() > 0) {
											_variables_B_a.put(vId, vSet_B);
											iter = vSet_B.vertices.iterator();
											while (iter.hasNext())
												_vertex_B_a[iter.next()].add(vId);
										}
									}
								}
							}
							if (rmp.variables_C != null) {
								for (int w = 0; w < n; w++)
									_vertex_C_a[w] = new ArrayList<Integer>();
								iterVarSet = rmp.variables_C.entrySet().iterator();
								while (iterVarSet.hasNext()) {
									HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
											.next();
									vSet = (VarSet) pairE.getValue();
									vId = pairE.getKey();
									if (vSet != null && vSet.vertices.size() > 0) {
										vSet_C = new VarSet();
										iter = vSet.vertices.iterator();
										while (iter.hasNext()) {
											v = iter.next();
											if (_uncovered_a[v])
												vSet_C.vertices.add(v);
										}
										vSet_C.size = vSet_C.vertices.size();
										if (vSet_C != null && vSet_C.vertices.size() > 0
												&& vSet_C.vertices.size() <= best_nr_C) {
											_variables_C_a.put(vId, vSet_C);
											iter = vSet_C.vertices.iterator();
											while (iter.hasNext())
												_vertex_C_a[iter.next()].add(vId);
										}
									}
								}
							}
							rmp_a = new RMP(_uncovered_a, _covered_A0_a, _covered_B0_a, _covered_C0_a, _eligible_A_a,
									_eligible_B_a, _nr_A0_a, _nr_B0_a, _nr_C0_a, _nr_eligible_A_a, _nr_eligible_B_a,
									_variables_A_a, _variables_B_a, _variables_C_a, _vertex_A_a, _vertex_B_a,
									_vertex_C_a);
							stackRMP.push(rmp_a);

						} else {// _nr_eligible_A_a == 0 or _nr_eligible_B_a == 0
							if (_nr_eligible_A_a == 0 && _nr_eligible_B_a > 0) {// all the vertices in V' are B
								// eligible;
								if (_nr_C0_a < best_nr_C) {
									newIntSolution = new int[n];
									for (int w = 0; w < n; w++) {
										if (_covered_A0_a[w])
											newIntSolution[w] = 1;
										if (_covered_B0_a[w] || _eligible_B_a[w])
											newIntSolution[w] = 2;
										if (_covered_C0_a[w])
											newIntSolution[w] = 3;
									}
									if (checkIntegerSolution(newIntSolution))
										bestSolution = newIntSolution;
								}
							}
							if (_nr_eligible_A_a > 0 && _nr_eligible_B_a == 0) {// all the vertices in V' are A
																				// eligible;
								if (_nr_C0_a < best_nr_C) {
									newIntSolution = new int[n];
									for (int w = 0; w < n; w++) {
										if (_covered_A0_a[w] || _eligible_A_a[w])
											newIntSolution[w] = 1;
										if (_covered_B0_a[w])
											newIntSolution[w] = 2;
										if (_covered_C0_a[w])
											newIntSolution[w] = 3;
									}
									if (checkIntegerSolution(newIntSolution))
										bestSolution = newIntSolution;
								}
							}
							if (_nr_eligible_A_a == 0 && _nr_eligible_B_a == 0) {// there are no remaining eligible
																					// vertices;
								if (_nr_C0_a < best_nr_C) {
									newIntSolution = new int[n];
									for (int w = 0; w < n; w++) {
										if (_covered_A0_a[w])
											newIntSolution[w] = 1;
										if (_covered_B0_a[w])
											newIntSolution[w] = 2;
										if (_covered_C0_a[w])
											newIntSolution[w] = 3;
									}
									if (checkIntegerSolution(newIntSolution))
										bestSolution = newIntSolution;
								}
							}
						}
					}
				}
				// B branch:
				if (rmp.eligible_B[u] && rmp.nr_B0 <= b - 1) {
					eligible_b = true;
					boolean[] _uncovered_b = new boolean[n], _covered_A0_b = new boolean[n],
							_covered_B0_b = new boolean[n], _covered_C0_b = new boolean[n],
							_eligible_A_b = new boolean[n], _eligible_B_b = new boolean[n];

					for (int w = 0; w < n; w++) {
						_uncovered_b[w] = rmp.uncovered[w];
						_covered_A0_b[w] = rmp.covered_A0[w];
						_covered_B0_b[w] = rmp.covered_B0[w];
						_covered_C0_b[w] = rmp.covered_C0[w];
						_eligible_A_b[w] = rmp.eligible_A[w];
						_eligible_B_b[w] = rmp.eligible_B[w];
					}
					_uncovered_b[u] = false;
					_covered_B0_b[u] = true;
					_eligible_B_b[u] = false;

					_nr_A0_b = rmp.nr_A0;
					_nr_B0_b = rmp.nr_B0 + 1;
					_nr_C0_b = rmp.nr_C0;

					_nr_eligible_A_b = rmp.nr_eligible_A;
					_nr_eligible_B_b = rmp.nr_eligible_B - 1;
					if (_eligible_A_b[u])
						_nr_eligible_A_b--;
					_eligible_A_b[u] = false;

					for (int w = 0; w < n; w++)
						if (u_neighborhood[w] && _eligible_A_b[w]) {
							_eligible_A_b[w] = false;
							_nr_eligible_A_b--;
						}

					for (int w = 0; w < n; w++)
						if (!_eligible_A_b[w] && !_eligible_B_b[w] && _uncovered_b[w]) {
							_uncovered_b[w] = false;
							_covered_C0_b[w] = true;
							_nr_C0_b++;
						}
					if (_nr_C0_b < best_nr_C) {
						if (_nr_eligible_A_b > 0 && _nr_eligible_B_b > 0) {

							HashMap<Integer, VarSet> _variables_A_b = new HashMap<Integer, VarSet>(),
									_variables_B_b = new HashMap<Integer, VarSet>(),
									_variables_C_b = new HashMap<Integer, VarSet>();

							ArrayList<Integer>[] _vertex_A_b = new ArrayList[n], _vertex_B_b = new ArrayList[n],
									_vertex_C_b = new ArrayList[n];

							if (rmp.variables_A != null && _nr_eligible_A_b > 0) {
								for (int w = 0; w < n; w++)
									_vertex_A_b[w] = new ArrayList<Integer>();
								iterVarSet = rmp.variables_A.entrySet().iterator();
								while (iterVarSet.hasNext()) {
									HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
											.next();
									vSet = (VarSet) pairE.getValue();
									vId = pairE.getKey();
									if (vSet != null && vSet.vertices.size() >= lowerBound_AB) {
										vSet_A = new VarSet();
										iter = vSet.vertices.iterator();
										while (iter.hasNext()) {
											v = iter.next();
											if (_eligible_A_b[v])
												vSet_A.vertices.add(v);
										}
										vSet_A.size = vSet_A.vertices.size();
										if (vSet_A != null && vSet_A.vertices.size() >= lowerBound_AB
												&& vSet_A.vertices.size() > 0) {
											_variables_A_b.put(vId, vSet_A);
											iter = vSet_A.vertices.iterator();
											while (iter.hasNext())
												_vertex_A_b[iter.next()].add(vId);
										}
									}
								}
							}
							if (rmp.variables_B != null && _nr_eligible_B_b > 0) {
								for (int w = 0; w < n; w++)
									_vertex_B_b[w] = new ArrayList<Integer>();
								iterVarSet = rmp.variables_B.entrySet().iterator();
								while (iterVarSet.hasNext()) {
									HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
											.next();
									vSet = (VarSet) pairE.getValue();
									vId = pairE.getKey();
									if (vSet != null && vSet.vertices.size() >= lowerBound_AB) {
										vSet_B = new VarSet();
										iter = vSet.vertices.iterator();
										while (iter.hasNext()) {
											v = iter.next();
											if (_eligible_B_b[v])
												vSet_B.vertices.add(v);
										}
										vSet_B.size = vSet_B.vertices.size();
										if (vSet_B != null && vSet_B.vertices.size() >= lowerBound_AB
												&& vSet_B.vertices.size() > 0) {
											_variables_B_b.put(vId, vSet_B);
											iter = vSet_B.vertices.iterator();
											while (iter.hasNext())
												_vertex_B_b[iter.next()].add(vId);
										}
									}
								}
							}
							if (rmp.variables_C != null) {
								for (int w = 0; w < n; w++)
									_vertex_C_b[w] = new ArrayList<Integer>();
								iterVarSet = rmp.variables_C.entrySet().iterator();
								while (iterVarSet.hasNext()) {
									HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet
											.next();
									vSet = (VarSet) pairE.getValue();
									vId = pairE.getKey();
									if (vSet != null && vSet.vertices.size() > 0) {
										vSet_C = new VarSet();
										iter = vSet.vertices.iterator();
										while (iter.hasNext()) {
											v = iter.next();
											if (_uncovered_b[v])
												vSet_C.vertices.add(v);
										}
										vSet_C.size = vSet_C.vertices.size();
										if (vSet_C != null && vSet_C.vertices.size() > 0
												&& vSet_C.vertices.size() <= best_nr_C) {
											_variables_C_b.put(vId, vSet_C);
											iter = vSet_C.vertices.iterator();
											while (iter.hasNext())
												_vertex_C_b[iter.next()].add(vId);
										}
									}
								}
							}
							rmp_b = new RMP(_uncovered_b, _covered_A0_b, _covered_B0_b, _covered_C0_b, _eligible_A_b,
									_eligible_B_b, _nr_A0_b, _nr_B0_b, _nr_C0_b, _nr_eligible_A_b, _nr_eligible_B_b,
									_variables_A_b, _variables_B_b, _variables_C_b, _vertex_A_b, _vertex_B_b,
									_vertex_C_b);
							stackRMP.push(rmp_b);

						} else {// _nr_eligible_A_b == 0 or _nr_eligible_B_b == 0
							if (_nr_eligible_A_b == 0 && _nr_eligible_B_b > 0) {// all the vertices in V' are B
								// eligible;
								if (_nr_C0_b < best_nr_C) {
									newIntSolution = new int[n];
									for (int w = 0; w < n; w++) {
										if (_covered_A0_b[w])
											newIntSolution[w] = 1;
										if (_covered_B0_b[w] || _eligible_B_b[w])
											newIntSolution[w] = 2;
										if (_covered_C0_b[w])
											newIntSolution[w] = 3;
									}
									if (checkIntegerSolution(newIntSolution))
										bestSolution = newIntSolution;
								}
							}
							if (_nr_eligible_A_b > 0 && _nr_eligible_B_b == 0) {// all the vertices in V' are A
																				// eligible;
								if (_nr_C0_b < best_nr_C) {
									newIntSolution = new int[n];
									for (int w = 0; w < n; w++) {
										if (_covered_A0_b[w] || _eligible_A_b[w])
											newIntSolution[w] = 1;
										if (_covered_B0_b[w])
											newIntSolution[w] = 2;
										if (_covered_C0_b[w])
											newIntSolution[w] = 3;
									}
									if (checkIntegerSolution(newIntSolution))
										bestSolution = newIntSolution;
								}
							}
							if (_nr_eligible_A_b == 0 && _nr_eligible_B_b == 0) {// there are no remaining eligible
																					// vertices;
								if (_nr_C0_b < best_nr_C) {
									newIntSolution = new int[n];
									for (int w = 0; w < n; w++) {
										if (_covered_A0_b[w])
											newIntSolution[w] = 1;
										if (_covered_B0_b[w])
											newIntSolution[w] = 2;
										if (_covered_C0_b[w])
											newIntSolution[w] = 3;
									}
									if (checkIntegerSolution(newIntSolution))
										bestSolution = newIntSolution;
								}
							}
						}
					}
				}

			} else {// rmp.nr_A0 == b && rmp.nr_B0 == b
				if (rmp.nr_C0 < best_nr_C) {
					newIntSolution = new int[n];
					for (int w = 0; w < n; w++) {
						if (rmp.covered_A0[w])
							newIntSolution[w] = 1;
						if (rmp.covered_B0[w])
							newIntSolution[w] = 2;
						if (rmp.covered_C0[w])
							newIntSolution[w] = 3;
					}
					if (checkIntegerSolution(newIntSolution))
						bestSolution = newIntSolution;
				}
			}
		}
	}

	public static boolean checkIntegerSolution(int[] solution) {
		int card_a = 0, card_b = 0, card_c;
		for (int u = 0; u < n; u++) {
			if (solution[u] == 1)
				card_a++;
			if (solution[u] == 2)
				card_b++;
			for (int v = 0; v < n; v++)
				if (solution[u] == 1 && solution[v] == 2 && adjacency[u][v])
					return false;
		}
		if (card_a > b || card_a == 0 || card_b > b || card_b == 0)
			return false;
		correctSolution = true;
		card_c = n - card_a - card_b;
		if (card_c < best_nr_C && card_a > 0 && card_b > 0) {
			best_nr_C = card_c;
			for (int u = 0; u < n; u++)
				bestSolution[u] = solution[u];
		}
		return true;
	}

	public static int[] rmpIntegerSolutionFound(GRBModel model, RMP rmp) throws GRBException {
		// The current rmp was solved and an integer solution was found.
		Iterator<Entry<Integer, VarSet>> iterVarSet;
		Iterator<Integer> iterVertex;
		double varValue;
		VarSet varSet;
		int nrC0 = 0, nrA0 = 0, nrB0 = 0, v;
		int[] newSolution = new int[n];
		for (int u = 0; u < n; u++) {
			if (rmp.covered_A0[u]) {
				newSolution[u] = 1;
				nrA0++;
			}
			if (rmp.covered_B0[u]) {
				newSolution[u] = 2;
				nrB0++;
			}
			if (rmp.covered_C0[u]) {
				newSolution[u] = 3;
				nrC0++;
			}
		}

		if (rmp.variables_A != null) {
			iterVarSet = rmp.variables_A.entrySet().iterator();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				varValue = model.getVarByName("x_" + pairE.getKey()).get(GRB.DoubleAttr.X);
				if (varValue >= 1 - precisionVar) {
					varSet = pairE.getValue();
					iterVertex = varSet.vertices.iterator();
					while (iterVertex.hasNext()) {
						v = iterVertex.next();
						nrA0++;
						if (newSolution[v] != 0) {
							// System.out.println("Integer solution conflict on A: " + newSolution[v] + ", "
							// + v);

						} else
							newSolution[v] = 1;
					}
				}
			}
		}

		if (rmp.variables_B != null) {
			iterVarSet = rmp.variables_B.entrySet().iterator();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				varValue = model.getVarByName("y_" + pairE.getKey()).get(GRB.DoubleAttr.X);
				if (varValue >= 1 - precisionVar) {
					varSet = pairE.getValue();
					iterVertex = varSet.vertices.iterator();
					while (iterVertex.hasNext()) {
						v = iterVertex.next();
						nrB0++;
						if (newSolution[v] != 0) {
							// System.out.println("Integer solution conflict on B: " + newSolution[v] + ", "
							// + v);
						} else
							newSolution[v] = 2;
					}
				}
			}
		}

		if (rmp.variables_C != null) {
			iterVarSet = rmp.variables_C.entrySet().iterator();
			while (iterVarSet.hasNext()) {
				HashMap.Entry<Integer, VarSet> pairE = (HashMap.Entry<Integer, VarSet>) iterVarSet.next();
				varValue = model.getVarByName("z_" + pairE.getKey()).get(GRB.DoubleAttr.X);
				if (varValue >= 1 - precisionVar) {
					varSet = pairE.getValue();
					iterVertex = varSet.vertices.iterator();
					while (iterVertex.hasNext()) {
						v = iterVertex.next();
						nrC0++;
						if (newSolution[v] != 0) {
							// System.out.println("Integer solution conflict on C: " + newSolution[v] + ", "
							// + v);
						} else
							newSolution[v] = 3;
					}
				}
			}
		}

		checkIntegerSolution(newSolution);
		if (nrC0 < best_nr_C && nrA0 > 0 && nrB0 > 0) {
			best_nr_C = nrC0;
			System.out.println("Better solution!" + " nr_a = " + nrA0 + ", nr_b = " + nrB0 + ", nr_c = " + nrC0);
			return newSolution;
		}
		return null;
	}

	public static RMP buildInitialRMP() {
		int a = 0, c = 0, h, min = n, u_0 = -1, card_a = 0, card_b = 0, card;
		boolean[] characteristic_A = new boolean[n], characteristic_B = new boolean[n],
				characteristic_C = new boolean[n], _uncovered = new boolean[n], _covered_A0 = new boolean[n],
				_covered_B0 = new boolean[n], _covered_C0 = new boolean[n], _eligible_A = new boolean[n],
				_eligible_B = new boolean[n];
		HashMap<Integer, VarSet> _variables_A = new HashMap<Integer, VarSet>(),
				_variables_B = new HashMap<Integer, VarSet>(), _variables_C = new HashMap<Integer, VarSet>();
		ArrayList<Integer>[] _vertex_A = new ArrayList[n], _vertex_B = new ArrayList[n], _vertex_C = new ArrayList[n];
		for (int u = 0; u < n; u++) {
			_vertex_A[u] = new ArrayList<Integer>();
			_vertex_B[u] = new ArrayList<Integer>();
			_vertex_C[u] = new ArrayList<Integer>();
		}

		for (int u = 0; u < n; u++) {
			_uncovered[u] = true;
			_eligible_A[u] = true;
			_eligible_B[u] = true;
			if (degrees[u] < n - 1)
				u_0 = u;
		}

		characteristic_A[u_0] = true;
		a++;
		for (int v = 0; v < n; v++)
			if (adjacency[u_0][v]) {
				characteristic_C[v] = true;
				c++;
			}
		while (a + c < n - b) {
			min = n;
			for (int u = 0; u < n; u++) {
				if (!characteristic_A[u]) {
					h = 0;
					for (int v = 0; v < n; v++)
						if (adjacency[u][v] && !characteristic_A[v] && !characteristic_C[v])
							h++;
					if (h < min) {
						min = h;
						u_0 = u;
					}
				}
			}
			characteristic_A[u_0] = true;
			a++;
			c = 0;
			for (int v = 0; v < n; v++)
				characteristic_C[v] = false;
			for (int v = 0; v < n; v++)
				if (!characteristic_A[v])
					for (int u = 0; u < n; u++)
						if (characteristic_A[u] && adjacency[u][v]) {
							characteristic_C[v] = true;
							c++;
							break;
						}
		}

		for (int v = 0; v < n; v++) {
			if (characteristic_A[v])
				card_a++;
			else if (!characteristic_C[v]) {
				characteristic_B[v] = true;
				card_b++;
			}
		}
		// |A| >= |B| for asymmetry reasons
		if (card_a < card_b) {
			boolean[] aux = new boolean[n];
			for (int v = 0; v < n; v++) {
				aux[v] = characteristic_B[v];
				characteristic_B[v] = characteristic_A[v];
				characteristic_A[v] = aux[v];
			}
			card = card_b;
			card_b = card_a;
			card_a = card;
		}
		VarSet vSet_A = new VarSet(), vSet_B = new VarSet(), vSet_C = new VarSet();
		for (int v = 0; v < n; v++) {
			if (characteristic_A[v]) {
				vSet_A.vertices.add(v);
				_vertex_A[v].add(0);
				bestSolution[v] = 1;
			} else if (characteristic_C[v]) {
				vSet_C.vertices.add(v);
				_vertex_C[v].add(2);
				bestSolution[v] = 3;
			} else {
				vSet_B.vertices.add(v);
				_vertex_B[v].add(1);
				bestSolution[v] = 2;
			}
		}
		vSet_A.size = vSet_A.vertices.size();
		vSet_B.size = vSet_B.vertices.size();
		vSet_C.size = vSet_C.vertices.size();
		_variables_A.put(0, vSet_A);
		_variables_B.put(1, vSet_B);
		_variables_C.put(2, vSet_C);
		varId = 3;
		best_nr_C = c;
		System.out.println("Initially |C| = " + best_nr_C);
		init_nr_C = best_nr_C;
		vSet_A.toString();
		vSet_B.toString();
		vSet_C.toString();
		checkSolution2(bestSolution);

		return new RMP(_uncovered, _covered_A0, _covered_B0, _covered_C0, _eligible_A, _eligible_B, 0, 0, 0, n, n,
				_variables_A, _variables_B, _variables_C, _vertex_A, _vertex_B, _vertex_C);
	}

	public static boolean checkSolution1(boolean[][] solution) {
		int nrA = 0, nrB = 0;
		for (int u = 0; u < n; u++)
			if (solution[0][u])
				for (int v = 0; v < n; v++)
					if (solution[1][v] && adjacency[u][v]) {
						System.out.println("Conflict, u = " + u + ", v = " + v);
						return false;
					}
		for (int u = 0; u < n; u++) {
			if (solution[0][u])
				nrA++;
			if (solution[1][u])
				nrB++;
		}
		if (nrA > b || nrA == 0 || nrB > b || nrB == 0)
			return false;
		return true;
	}

	public static void checkSolution2(int[] solution) {
		int nrA = 0, nrB = 0;
		for (int u = 0; u < n; u++)
			if (solution[u] == 1)
				for (int v = 0; v < n; v++)
					if (solution[v] == 2 && adjacency[u][v])
						System.out.println("Conflict, u = " + u + ", v = " + v);
		for (int u = 0; u < n; u++)
			if (solution[u] == 1)
				nrA++;
		for (int u = 0; u < n; u++)
			if (solution[u] == 2)
				nrB++;

		if (nrA > b || nrA == 0)
			System.out.println("(A) Size conflict!");
		else if (nrB > b || nrB == 0)
			System.out.println("(B) Size conflict!");
		else
			System.out.println("Initial solution is a correct one.");
	}

	public static void solveRMP_Pool(RMP rmp, int noPb, double elapsed) throws GRBException {
		// solve the RMP and branch when appropriate
		double[] duals_alpha = new double[n], a_weight = new double[n], b_weight = new double[n],
				c_weight = new double[n];
		double dual_beta_a, dual_beta_b, dual_delta_a, dual_delta_b, dual_delta_c, dual_epsilon;
		long rmp_solving_time = System.nanoTime(), rmp_total_time = System.nanoTime();

		GRBModel _model = buildLPModel(rmp);
		rmp_solving_time = System.nanoTime();
		GRBConstr _constraint;
		int status = -1, vertex, step = 0, nr_A0 = rmp.nr_A0, nr_B0 = rmp.nr_B0, nr_C0 = rmp.nr_C0;
		boolean found = true, found_a, found_b, found_c;
		VarSet[] varSet_AA, varSet_BB, varSet_CC;
		while (found) {
			if (elapsed + (System.nanoTime() - rmp_total_time) * 1e-9 > time_max)
				return;
			step++;
			_model.set(GRB.IntParam.Method, 3);
			_model.set(GRB.IntParam.ConcurrentMethod, 0);
			_model.optimize();
			status = _model.get(GRB.IntAttr.Status);
			if (status == GRB.Status.OPTIMAL) {
				System.out.println("%%%%%% Rmp value = " + _model.get(DoubleAttr.ObjVal) + ", step " + step + "%%%%%%");
			} else
				System.out.println("%%%%%% No optimal solution %%%%%%");
			if (status == GRB.Status.OPTIMAL) {
				dual_beta_a = _model.getConstrByName("constraint_A").get(GRB.DoubleAttr.Pi);
				dual_beta_b = _model.getConstrByName("constraint_B").get(GRB.DoubleAttr.Pi);

				_constraint = _model.getConstrByName("one_subset_A");
				dual_delta_a = 0;
				if (_constraint != null) {
					dual_delta_a = _constraint.get(GRB.DoubleAttr.Pi);
				} else {
					System.out.println("dual_delta_a NULL!");
				}

				_constraint = _model.getConstrByName("one_subset_B");
				dual_delta_b = 0;
				if (_constraint != null) {
					dual_delta_b = _constraint.get(GRB.DoubleAttr.Pi);
				} else {
					System.out.println("dual_delta_b NULL!");
				}

				_constraint = _model.getConstrByName("one_subset_C");
				dual_delta_c = 0;
				if (_constraint != null) {
					dual_delta_c = _constraint.get(GRB.DoubleAttr.Pi);
				} else {
					System.out.println("dual_delta_c NULL!");
				}

				_constraint = _model.getConstrByName("asymmetric");
				dual_epsilon = 0;
				if (_constraint != null)
					dual_epsilon = _model.getConstrByName("asymmetric").get(GRB.DoubleAttr.Pi);

				found_a = found_b = found_c = true;

				for (int u = 0; u < n; u++)
					if (rmp.uncovered[u]) {
						duals_alpha[u] = _model.getConstrByName("vertex_" + u).get(GRB.DoubleAttr.Pi);
						c_weight[u] = duals_alpha[u] - 1;
					}

				for (int u = 0; u < n; u++) {
					if (rmp.eligible_A[u])
						a_weight[u] = duals_alpha[u] + dual_beta_a + dual_epsilon;
					if (rmp.eligible_B[u])
						b_weight[u] = duals_alpha[u] + dual_beta_b - dual_epsilon;
				}

				for (int u = 0; u < n; u++)
					if (rmp.eligible_A[u]) {
						for (int v = u + 1; v < n; v++)
							if (adjacency[u][v] && rmp.eligible_B[v]) {
								_constraint = _model.getConstrByName("A_B_edge_" + u + "_" + v);
								if (_constraint != null) {
									a_weight[u] += _constraint.get(GRB.DoubleAttr.Pi);// gamma_{uv}
									b_weight[v] += _constraint.get(GRB.DoubleAttr.Pi);// gamma_{uv}}
								}
							}
						for (int v = 0; v < u; v++)
							if (adjacency[u][v] && rmp.eligible_B[v]) {
								_constraint = _model.getConstrByName("B_A_edge_" + v + "_" + u);
								if (_constraint != null) {
									a_weight[u] += _constraint.get(GRB.DoubleAttr.Pi);// gamma_{vu}
									b_weight[v] += _constraint.get(GRB.DoubleAttr.Pi);// gamma_{vu}}
								}
							}
					}

				found_a = false;
				varSet_AA = solveSubProblem(a_weight, rmp.eligible_A, "A", nr_A0, -dual_delta_a);
				if (varSet_AA != null)
					for (int j = 0; j < varSet_AA.length; j++)
						if (varSet_AA[j] != null && varSet_AA[j].vertices.size() > 0) {
							addVarSet(_model, rmp, varSet_AA[j], "A");
							found_a = true;
						}

				found_b = false;
				varSet_BB = solveSubProblem(b_weight, rmp.eligible_B, "B", nr_B0, -dual_delta_b);
				if (varSet_BB != null)
					for (int j = 0; j < varSet_BB.length; j++)
						if (varSet_BB[j] != null && varSet_BB[j].vertices.size() > 0) {
							addVarSet(_model, rmp, varSet_BB[j], "B");
							found_b = true;
						}

				found_c = false;
				varSet_CC = solveSubProblem(c_weight, rmp.uncovered, "C", nr_C0, -dual_delta_c);
				if (varSet_CC != null)
					for (int j = 0; j < varSet_CC.length; j++)
						if (varSet_CC[j] != null && varSet_CC[j].vertices.size() > 0) {
							addVarSet(_model, rmp, varSet_CC[j], "C");
							found_c = true;
						}
				if (found_a == false && found_b == false && found_c == false)
					found = false;
			} else {
				if (status == GRB.Status.INFEASIBLE)
					System.out.println("RMPImproved infeasible!");
				if (status == GRB.Status.UNBOUNDED)
					System.out.println("RMPImproved unbounded!");
				if (status == GRB.Status.INF_OR_UNBD)
					System.out.println("RMPImproved infeasible or unbounded!");
				found = false;
			}
		}

		System.out.println("** Solving time for rmp: " + (System.nanoTime() - rmp_solving_time) * 1e-9 + " **");
		if (status == GRB.Status.OPTIMAL && _model.get(DoubleAttr.ObjVal) < best_nr_C) {
			vertex = branchOnVertex(_model, rmp);
			if (vertex != -1)
				branchNImproved(rmp, vertex);
			else {
				readSolution(rmp, _model);
				rmpIntegerSolutionFound(_model, rmp);
			}
		}
		if (noPb == 1) {
			rootTime = (System.nanoTime() - rmp_total_time) * 1e-9;
			rootNoOfVars = 0;
			if (rmp.variables_A != null && rmp.variables_A.size() > 0)
				rootNoOfVars = rmp.variables_A.size();
			if (rmp.variables_B != null && rmp.variables_B.size() > 0)
				rootNoOfVars += rmp.variables_B.size();
			if (rmp.variables_C != null && rmp.variables_C.size() > 0)
				rootNoOfVars += rmp.variables_C.size();
		}

		totalNoOfVars = 0;
		if (rmp.variables_A != null && rmp.variables_A.size() > 0)
			totalNoOfVars = rmp.variables_A.size();
		if (rmp.variables_B != null && rmp.variables_B.size() > 0)
			totalNoOfVars += rmp.variables_B.size();
		if (rmp.variables_C != null && rmp.variables_C.size() > 0)
			totalNoOfVars += rmp.variables_C.size();

		System.out.println("BEST c: " + best_nr_C);
		return;
	}

	public static LinkedList<Entry<Integer, Double>> weightsSort(double[] weight, final boolean order) {
		LinkedList<Entry<Integer, Double>> list = new LinkedList<>();
		for (int i = 0; i < n; i++)
			list.add(Map.entry(i, weight[i]));

		// Sorting the list based on values
		list.sort((o1, o2) -> order
				? o1.getValue().compareTo(o2.getValue()) == 0 ? o1.getKey().compareTo(o2.getKey())
						: o1.getValue().compareTo(o2.getValue())
				: o2.getValue().compareTo(o1.getValue()) == 0 ? o2.getKey().compareTo(o1.getKey())
						: o2.getValue().compareTo(o1.getValue()));
		return list;
	}

	public static VarSet[] solveSubProblem(double[] weight, boolean[] eligible, String type, int nr_A0_B0_C0,
			double minus_dual_delta) {
		if (Math.abs(minus_dual_delta) == 0)
			minus_dual_delta = 0;
		VarSet[] sets;
		VarSet varS;
		double setWeight = 0, weight_v;
		int v = 0, size = 0, index_first = -1, l, lowerBound_AB = n - b - best_nr_C + 1, index_last_zero = -1;
		LinkedList<Integer> addedVertex = new LinkedList<Integer>();
		LinkedList<Entry<Integer, Double>> list = weightsSort(weight, false);

		if (type != "C") {// type == "A" or "B"
			size = 0;
			index_first = -1;
			for (Entry<Integer, Double> entry : list) {
				weight_v = entry.getValue();
				size++;
				if (weight_v == 0)
					index_last_zero = size - 1;
				if (weight_v < 0)
					break;
			}

			size = 0;
			for (Entry<Integer, Double> entry : list) {
				weight_v = entry.getValue();
				v = entry.getKey();
				if (eligible[v] && size < b - nr_A0_B0_C0) {
					if (weight_v < 0
							&& (index_first >= 0 || setWeight + weight_v <= minus_dual_delta + precisionSetWeight))
						break;
					setWeight += weight_v;
					addedVertex.add(v);
					size++;
					if (setWeight > minus_dual_delta + precisionSetWeight && index_first == -1)
						index_first = size - 1;
				}
			}
			if (size < lowerBound_AB)// the possible sets are too small
				size = -1;
		} else {// type == "C"
			size = 0;
			index_first = -1;
			for (Entry<Integer, Double> entry : list) {
				weight_v = entry.getValue();
				v = entry.getKey();
				if (eligible[v] && size < best_nr_C - 1 - nr_A0_B0_C0) {
					if (weight_v < 0
							&& (index_first >= 0 || setWeight + weight_v <= minus_dual_delta + precisionSetWeight))
						break;
					setWeight += weight_v;
					addedVertex.add(v);
					size++;
					if (setWeight > minus_dual_delta + precisionSetWeight && size >= alpha_star - nr_A0_B0_C0
							&& index_first == -1)
						index_first = size - 1;
				}
				if (size >= best_nr_C - 1 - nr_A0_B0_C0)
					break;
			}
		}

		if (type != "C" && setWeight > minus_dual_delta + precisionSetWeight && index_first >= 0) {
			int local_size = 0;
			setWeight = 0;
			sets = new VarSet[1];
			varS = new VarSet();
			Iterator<Integer> iterA = addedVertex.iterator();
			while (iterA.hasNext()) {
				v = iterA.next();
				weight_v = weight[v];
				if (weight_v >= 0) {
					local_size++;
					setWeight += weight_v;
					varS.vertices.add(v);
					if (weight_v == 0 && local_size >= lowerBound_AB)
						break;
				}
			}
			varS.size = varS.vertices.size();
			sets[0] = varS;
			return sets;
		}

		if (type == "C" && index_first >= 0 && setWeight > minus_dual_delta + precisionSetWeight && size > 0) {
			sets = new VarSet[size - index_first];
			for (int h = 0; h < size - index_first; h++) {
				varS = new VarSet();
				Iterator<Integer> iterA = addedVertex.iterator();
				l = 0;
				while (iterA.hasNext()) {
					v = iterA.next();
					if (l <= index_first + h) {
						l++;
						varS.vertices.add(v);
					} else
						break;
				}
				varS.size = varS.vertices.size();
				sets[h] = varS;
			}
			return sets;
		}
		return null;
	}

	public static RMP addVarSet(GRBModel model, RMP rmp, VarSet set, String type) throws GRBException {
		// called only if set is not NULL;
		Iterator<Integer> iteratorV;
		int u, h = 0, q = 0;
		double[] coeffConstr;
		GRBConstr[] newConstr;
		GRBConstr _constraint;
		VarSet newSet = new VarSet();
		iteratorV = set.vertices.iterator();
		while (iteratorV.hasNext())
			newSet.vertices.add(iteratorV.next());
		newSet.size = newSet.vertices.size();

		if (type == "A") {
			_constraint = model.getConstrByName("constraint_A");
			if (_constraint != null)
				q++;
			_constraint = model.getConstrByName("one_subset_A");
			if (_constraint != null)
				q++;
			_constraint = model.getConstrByName("asymmetric");
			if (_constraint != null)
				q++;

			iteratorV = set.vertices.iterator();
			while (iteratorV.hasNext()) {
				u = iteratorV.next();
				for (int v = 0; v < n; v++)
					if (adjacency[u][v] && rmp.eligible_B[v])
						q++;
			}
			q += set.size;
			coeffConstr = new double[q];
			newConstr = new GRBConstr[q];

			rmp.variables_A.put(varId, newSet);
			iteratorV = set.vertices.iterator();
			while (iteratorV.hasNext()) {
				u = iteratorV.next();
				rmp.vertex_A[u].add(varId);
				newConstr[h] = model.getConstrByName("vertex_" + u);
				coeffConstr[h] = 1;
				h++;
				for (int v = 0; v < n; v++)
					if (adjacency[u][v] && rmp.eligible_B[v]) {
						if (u < v) {
							_constraint = model.getConstrByName("A_B_edge_" + u + "_" + v);
							if (_constraint != null) {
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							} else {
								_constraint = model.addConstr(new GRBLinExpr(), GRB.LESS_EQUAL, 1,
										"A_B_edge_" + u + "_" + v);
								model.update();
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							}
						}
						if (u > v) {
							_constraint = model.getConstrByName("B_A_edge_" + v + "_" + u);
							if (_constraint != null) {
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							} else {
								_constraint = model.addConstr(new GRBLinExpr(), GRB.LESS_EQUAL, 1,
										"B_A_edge_" + v + "_" + u);
								model.update();
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							}
						}
					}
			}
			_constraint = model.getConstrByName("constraint_A");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = set.size;
				h++;
			}
			_constraint = model.getConstrByName("one_subset_A");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = 1.0;
				h++;
			}
			_constraint = model.getConstrByName("asymmetric");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = set.size;
			}
			// System.out.println("Id (addVar, A): " + varId);
			model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, newConstr, coeffConstr, "x_" + varId);
			varId++;
			model.update();
		}
		if (type == "B") {
			_constraint = model.getConstrByName("constraint_B");
			if (_constraint != null)
				q++;
			_constraint = model.getConstrByName("one_subset_B");
			if (_constraint != null)
				q++;
			_constraint = model.getConstrByName("asymmetric");
			if (_constraint != null)
				q++;

			iteratorV = set.vertices.iterator();
			while (iteratorV.hasNext()) {
				u = iteratorV.next();
				for (int v = 0; v < n; v++)
					if (adjacency[u][v] && rmp.eligible_A[v])
						q++;
			}
			q += set.size;
			coeffConstr = new double[q];
			newConstr = new GRBConstr[q];

			rmp.variables_B.put(varId, newSet);
			iteratorV = set.vertices.iterator();
			while (iteratorV.hasNext()) {
				u = iteratorV.next();
				rmp.vertex_B[u].add(varId);
				newConstr[h] = model.getConstrByName("vertex_" + u);
				coeffConstr[h] = 1;
				h++;
				for (int v = 0; v < n; v++)
					if (adjacency[u][v] && rmp.eligible_A[v]) {
						coeffConstr[h] = 1;
						if (u < v) {
							_constraint = model.getConstrByName("B_A_edge_" + u + "_" + v);
							if (_constraint != null) {
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							} else {
								_constraint = model.addConstr(new GRBLinExpr(), GRB.LESS_EQUAL, 1,
										"B_A_edge_" + u + "_" + v);
								model.update();
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							}
						}
						if (u > v) {
							_constraint = model.getConstrByName("A_B_edge_" + v + "_" + u);
							if (_constraint != null) {
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							} else {
								_constraint = model.addConstr(new GRBLinExpr(), GRB.LESS_EQUAL, 1,
										"A_B_edge_" + v + "_" + u);
								model.update();
								newConstr[h] = _constraint;
								coeffConstr[h] = 1;
								h++;
							}
						}
					}
			}
			_constraint = model.getConstrByName("constraint_B");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = set.size;
				h++;
			}
			_constraint = model.getConstrByName("one_subset_B");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = 1.0;
				h++;
			}
			_constraint = model.getConstrByName("asymmetric");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = -set.size;
			}
			model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, newConstr, coeffConstr, "y_" + varId);
			varId++;
			model.update();
		}
		if (type == "C") {
			_constraint = model.getConstrByName("one_subset_C");
			if (_constraint != null)
				q++;
			q += set.size;
			coeffConstr = new double[q];
			newConstr = new GRBConstr[q];

			rmp.variables_C.put(varId, newSet);
			iteratorV = set.vertices.iterator();
			while (iteratorV.hasNext())
				while (iteratorV.hasNext()) {
					u = iteratorV.next();
					rmp.vertex_C[u].add(varId);
					newConstr[h] = model.getConstrByName("vertex_" + u);
					coeffConstr[h] = 1;
					h++;
				}
			_constraint = model.getConstrByName("one_subset_C");
			if (_constraint != null) {
				newConstr[h] = _constraint;
				coeffConstr[h] = 1.0;
			}
			model.addVar(0, GRB.INFINITY, set.size, GRB.CONTINUOUS, newConstr, coeffConstr, "z_" + varId);
			varId++;
			model.update();
		}
		return rmp;
	}

	public static void writeFinalData(int pbNo, double time) {
		createFile("../data/VSP/results/finalSolutions.txt");
		NumberFormat nrformat = NumberFormat.getInstance();
		nrformat.setMaximumFractionDigits(2);
		try (BufferedWriter writer = new BufferedWriter(
				new FileWriter("../data/VSP/results/finalSolutions.txt", true))) {
			writer.newLine();
			writer.newLine();
			String nLine = "Input file: " + fileName + ", " + n + " vertices" + ", " + m + " edges ("
					+ (n * (n - 1) / 2 - m) + " non-adjacencies)";
			writer.write(nLine);
			writer.newLine();

			writer.newLine();
			nLine = "******* Best separator size: " + bestSeparatorSize;
			System.out.println(nLine);
			writer.write(nLine);

			writer.newLine();
			nLine = "******* Initially separator size: " + init_nr_C;
			System.out.println(nLine);
			writer.write(nLine);

			writer.newLine();
			nLine = "******* Root objective: " + nrformat.format(rootObj);
			System.out.println(nLine);
			writer.write(nLine);
			System.out.println();
			writer.newLine();

			writer.newLine();
			nLine = "******* # of problems in the Branch&Price tree: " + pbNo;
			System.out.println(nLine);
			writer.write(nLine);
			System.out.println();
			writer.newLine();

			writer.newLine();
			nLine = "******* Solving time per problem: " + nrformat.format(time / pbNo) + "s";
			System.out.println(nLine);
			writer.write(nLine);

			writer.newLine();
			nLine = "******* Root solving time: " + nrformat.format(rootTime) + "s";
			System.out.println(nLine);
			writer.write(nLine);

			writer.newLine();
			nLine = "******* Total solving time: " + nrformat.format(time) + "s";
			System.out.println(nLine);
			writer.write(nLine);
			System.out.println();
			writer.newLine();

			writer.newLine();
			nLine = "******* Average # of variables: " + nrformat.format((double) totalNoOfVars / (double) pbNo);
			System.out.println(nLine);
			writer.write(nLine);

			writer.newLine();
			nLine = "******* Root # of variables: " + rootNoOfVars;
			System.out.println(nLine);
			writer.write(nLine);

			writer.newLine();
			nLine = "******* Total # of variables: " + totalNoOfVars;
			System.out.println(nLine);
			writer.write(nLine);

		} catch (IOException ex) {
			ex.printStackTrace();
		}
	}

	public static void problemInfo(RMP rmp) {
		System.out.println(
				"nr A_0: " + rmp.nr_A0 + ", nr B_0: " + rmp.nr_B0 + ", nr C_0: " + rmp.nr_C0 + ",best C: " + best_nr_C);
		System.out.print("A0: ");
		for (int u = 0; u < n; u++)
			if (rmp.covered_A0[u])
				System.out.print(u + ", ");
		System.out.println();
		System.out.print("B0: ");
		for (int u = 0; u < n; u++)
			if (rmp.covered_B0[u])
				System.out.print(u + ", ");
		System.out.println();
		System.out.print("C0: ");
		for (int u = 0; u < n; u++)
			if (rmp.covered_C0[u])
				System.out.print(u + ", ");
	}

	public static void writeBestSolution(int[] numbers, String filename) {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
			for (int num : numbers) {
				writer.write(Integer.toString(num));
				writer.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static long getSeparatorSize() {
		System.out.println("Best C (count): " + Arrays.stream(bestSolution).filter(x -> x == 3).count());
		System.out.println("BEST C (count): " + best_nr_C);
		return Arrays.stream(bestSolution).filter(x -> x == 3).count();
	}

	public static void branchAndPrice(int n_max) throws GRBException, IOException {
		int pbNo = 0;
		double elapsed = 0.0, time;
		long pb_time = System.nanoTime(), total_time = System.nanoTime();
		String path = "../data/PartitionCol/results/" + fileName;

		stackRMP = new Stack<RMP>();
		createDir(path);

		RMP rmp = buildInitialRMP();
		stackRMP.add(rmp);
		while (!stackRMP.isEmpty() && pbNo < n_max && elapsed < time_max) {
			pb_time = System.nanoTime();
			pbNo++;
			System.out.println();
			System.out.println("** Starting problem " + pbNo);
			rmp = stackRMP.pop();
			elapsed = (System.nanoTime() - total_time) * 1e-9;
			solveRMP_Pool(rmp, pbNo, elapsed);
			System.out.println(
					"** Solving time for problem " + pbNo + ": " + (System.nanoTime() - pb_time) * 1e-9 + " **");
			elapsed = (System.nanoTime() - total_time) * 1e-9;
			System.out.println("******* Total time: " + elapsed);
			System.out.println("** Finishing problem " + pbNo);

			long currentBestSeparator = getSeparatorSize();
			if (currentBestSeparator < bestSeparatorSize) {
				bestSeparatorSize = (int) currentBestSeparator;
			}
			System.out.println("Best C: " + bestSeparatorSize);
		}
		time = (System.nanoTime() - total_time) * 1e-9;
		writeFinalData(pbNo, time);
	}

	// (1) pop the problem from the stack;
	// (2) build the model;
	// (3) solve the model;
	// (4) branch if necessary.

	public static void main(String[] args) throws GRBException, IOException, InterruptedException {

		int n_max = 1_000_000;
		fileNameInstances = "instances1.txt";

		ArrayList<String> listInst = createListOfInstances(fileNameInstances);

		Iterator<String> iteratorInstance = listInst.iterator();
		while (iteratorInstance.hasNext()) {
			fileName = iteratorInstance.next();
			// readGraphFromFile(fileName);
			readGraphFromFileTable(fileName);
			b = (int) Math.floor((double) (2 * n) / 3);
			totalNoOfVars = 0;
			rootNoOfVars = 0;
			best_nr_C = n;
			bestSolution = new int[n];
			bestSeparatorSize = n;
			init_nr_C = 0;

			createDir(path_results + fileName);
			branchAndPrice(n_max);
			writeBestSolution(bestSolution, path_solutions + fileName);

			System.out.println();
			System.out.println("** END PROBLEM **");
			System.out.println();
		}
	}
}