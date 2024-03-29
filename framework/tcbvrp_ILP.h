#ifndef __TCBVRP_ILP__H__
#define __TCBVRP_ILP__H__

#include "Tools.h"
#include "Instance.h"
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN

class tcbvrp_ILP
{

	typedef IloArray<IloIntVarArray> IloIntVarArray2;
	typedef IloArray<IloBoolVarArray> IloBoolVarArray2;

	typedef IloArray<IloBoolVarArray>	BoolVarMatrix;
	typedef IloArray<BoolVarMatrix>		BoolVar3Matrix;		//3D Bool array
	typedef IloArray<BoolVar3Matrix>	BoolVar4Matrix;		//4D Bool array

	typedef IloArray<IloNumVarArray>	NumVarMatrix;
	typedef IloArray<NumVarMatrix>		NumVar3Matrix;		//3D Numerical array
	typedef IloArray<NumVar3Matrix>		NumVar4Matrix;		//4D Numerical array

private:

	Instance& instance;
	string model_type;

	unsigned int n; // Number of Stations + Depot
	unsigned int a; // Number of arcs
	unsigned int m; // Number of Vehicles
	unsigned int T; // Time budget

	IloEnv env;
	IloModel model;
	IloCplex cplex;

	void initCPLEX();
	void setCPLEXParameters();

	void initConstraints(BoolVar3Matrix var_t,IloBoolVarArray var_r);
	void initDecisionVars(BoolVar3Matrix &var_t, IloBoolVarArray &var_r);
	void initObjectiveFunction(BoolVar3Matrix var_t);
	void modelSCF();
	void modelMCF();
	void modelMTZ();

public:

	tcbvrp_ILP( Instance& _instance, string _model_type);
	~tcbvrp_ILP();
	void solve();

};

#endif //__TCBVRP_ILP__H__
