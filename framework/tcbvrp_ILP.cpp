#include "tcbvrp_ILP.h"

tcbvrp_ILP::tcbvrp_ILP( Instance& _instance, string _model_type) :
instance( _instance ), model_type( _model_type )
{
	//Number of stations + depot
	n = instance.n;
	//Number of edges
	a = instance.nArcs;
	//Max. Number of vehicles
	m = instance.m;
	//Max. time limit
	T = instance.T;
}

tcbvrp_ILP::~tcbvrp_ILP()
{
	// free CPLEX resources
	cplex.end();
	model.end();
	env.end();
}

void tcbvrp_ILP::solve()
{
	try {
		// initialize CPLEX solver
		env = IloEnv();
		model = IloModel( env );

		// add model-specific constraints
		if( model_type == "scf" )
			modelSCF();
		else if( model_type == "mcf" )
			modelMCF();
		else if( model_type == "mtz" )
			modelMTZ();

		// build model
		cplex = IloCplex( model );

		// export model to a text file
		//cplex.exportModel( "model.lp" );

		// set parameters
		setCPLEXParameters();

		// solve model
		cout << "Calling CPLEX solve ...\n";
		cplex.solve();
		cout << "CPLEX finished." << "\n\n";
		cout << "CPLEX status: " << cplex.getStatus() << "\n";
		cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << "\n";
		cout << "Objective value: " << cplex.getObjValue() << "\n";
		cout << "CPU time: " << Tools::CPUtime() << "\n\n";

		//Copied from: http://orinanobworld.blogspot.co.at/2010/08/iterating-over-cplex-model-sequel.html
		for (IloIterator<IloNumVar>  it(env); it.ok(); ++it) {
			try {
			IloNumVar ext = *it;
			if(cplex.getValue(ext)!=0)
				cout << ext.getName() << ": " << cplex.getValue(ext) << "\n";
			} catch (IloException &e){
				//cerr << "Exception for variable " << e << "\n";
			}
		}
	}
	catch( IloException& e ) {
		cerr << "tcbvrp_ILP: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cerr << "tcbvrp_ILP: unknown exception.\n";
		exit( -1 );
	}
}

// ----- private methods -----------------------------------------------

void tcbvrp_ILP::setCPLEXParameters()
{
	// print every line of node-log and give more details
	cplex.setParam( IloCplex::MIPInterval, 1 );
	cplex.setParam( IloCplex::MIPDisplay, 2 );

	// deactivate CPLEX general-purpose cuts
	//	cplex.setParam( IloCplex::EachCutLim, 0 );
	//	cplex.setParam( IloCplex::FracCuts, -1 );

	// only use a single thread
	cplex.setParam( IloCplex::Threads, 1 );
	// set time limit for cplex (in seconds)
	cplex.setParam( IloCplex::TiLim, 3600);
}

void tcbvrp_ILP::initObjectiveFunction(BoolVar3Matrix var_t)
{
	IloExpr objFunction(env);
	for(int i=0;i<instance.m;i++)
	{
		for(int j=0; j< instance.n; j++)
		{
			for(int k=0; k< instance.n; k++)
			{
				objFunction += var_t[i][j][k] * instance.getDistance(j, k);
			}
		}
	}

	model.add(IloMinimize(env, objFunction));
	objFunction.end();
}

void tcbvrp_ILP::initDecisionVars(BoolVar3Matrix &var_t, IloBoolVarArray &var_r)
{
	/*
	 * t(i,j,k) is 1 if the arc from (j,k) is used by the tour i
	 */

	for(int i=0; i < instance.m; i++)
	{
		var_t[i] = BoolVarMatrix(env, instance.n);
		for(int j=0; j < instance.n; j++)
		{
			var_t[i][j] = IloBoolVarArray(env, instance.n);
			for(int k=0; k < instance.n; k++)
			{
				var_t[i][j][k] = IloBoolVar(env, Tools::indicesToString( "t_", i, j, k).c_str());
			}
		}
	}

	/*
	 * additional (continuous) variables n(i) represent the number of nodes in the tour i
	 */

	for(int i=0; i < instance.m; i++)
	{
		var_r[i] = IloBoolVar(env, Tools::indicesToString( "r_", i).c_str());
	}
}

void tcbvrp_ILP::initConstraints(BoolVar3Matrix var_t,IloBoolVarArray var_r){

	/*
	* var_f is true if there is an outgoing route from the originator
	*/

	 for(int i=0;i<instance.m;i++)
	 {
	 	IloExpr exprExpr(env);
	 	for(int k=1; k < instance.n; k++)
	 	{
	 		exprExpr += var_t[i][0][k];
	 	}
	 	model.add(exprExpr == var_r[i]);
	 	exprExpr.end();
	 }

	 /*
	 * there are only other arcs if there is an outgoing arc from the originator
	 */

	 for(int i=0;i<instance.m;i++)
	 {
	 	for(int j=1; j < instance.n; j++)
	 	{
	 		for(int k=1; k < instance.n; k++)
	 		{
	 			model.add(var_t[i][j][k] <= var_r[i]);
	 		}
	 	}
	 }

	/*
	 * a supply node is not allowed to got to a supply node
	 * a demand node is not allowed to got to a demand node
	 * the originator is not allowed to got to a demand node
	 * a supply node is not allowed to got to the originator
	 * no self loops are allowed
	 */

	 for(int k=0; k < instance.m; k++)
	 {
	 	for(int i=0;i<instance.n;i++)
	 	{
	 		for(int j=0; j< instance.n; j++)
	 		{
	 			if(instance.isSupplyNode(i) && instance.isSupplyNode(j))
	 			{
	 				model.add(var_t[k][i][j] == 0);
	 			}
	 			if(instance.isDemandNode(i) && instance.isDemandNode(j))
	 			{
	 				model.add(var_t[k][i][j] == 0);
	 			}
	 			if(i==0 && instance.isDemandNode(j))
	 			{
	 				model.add(var_t[k][i][j] == 0);
	 			}
	 			if(instance.isSupplyNode(i) && j==0)
	 			{
	 				model.add(var_t[k][i][j] == 0);
	 			}
	 			if(i==j)
	 			{
	 				model.add(var_t[k][i][j] == 0);
	 			}
	 		}
	 	}
	 }

	/*
	 * Each demand node has to have an outgoing arc which goes to a supply node or the originator
	 */

	for(int j=1; j< instance.n; j++)
	{
		if(instance.isDemandNode(j))
		{
			IloExpr toSupplyExpr(env);
			for(int i=0;i<instance.m;i++)
			{
				for(int k=0; k < instance.n; k++)
				{
					if(instance.isSupplyNode(k) || k == 0)
					{
						toSupplyExpr += var_t[i][j][k];
					}
				}
			}
			model.add(toSupplyExpr == 1);
			toSupplyExpr.end();
		}
	}

	/*
	 * Each supply node can only go to at most one demand node
	 */

	for(int j=1; j< instance.n; j++)
	{
		if(instance.isSupplyNode(j))
		{
			IloExpr myExpr1(env);
			for(int i=0;i<instance.m;i++)
			{
				for(int k=0; k < instance.n; k++)
				{
					if(instance.isDemandNode(k))
					{
						myExpr1 += var_t[i][j][k];
					}
				}
			}
			model.add(myExpr1 <= 1);
			myExpr1.end();
		}
	}

	/*
	 * The originator is not allowed to have more than m outgoing arcs
	 */

	IloExpr myExpr3(env);
	for(int i=0;i<instance.m;i++)
	{
		for(int k=0; k < instance.n; k++)
		{
			myExpr3 += var_t[i][0][k];
		}
	}
	model.add(myExpr3 <= instance.m);
	myExpr3.end();

	/*
	 * If the ingoing arc in one node is from a tour the outgoing arc has to be from the same tour
	 */

	for(int j=0; j< instance.n; j++)
	{
		for(int i=0;i<instance.m;i++)
		{
			IloExpr myExpr4(env);
			IloExpr myExpr5(env);
			for(int k=0; k < instance.n; k++)
			{
				myExpr4 += var_t[i][k][j];
				myExpr5 += var_t[i][j][k];
			}
			model.add(myExpr4 == myExpr5);
			myExpr4.end();
			myExpr5.end();
		}
	}

	/*
	 * A tour must be finished under the maximum time
	 */

	for(int i=0;i<instance.m;i++)
	{
		IloExpr maxTimeExpr(env);
		for(int j=0;j<instance.n;j++)
		{
			for(int k=0; k < instance.n; k++)
			{
				maxTimeExpr += var_t[i][j][k] * instance.getDistance(j, k);
			}
		}
		model.add(maxTimeExpr <= instance.T);
		maxTimeExpr.end();
	}
}

void tcbvrp_ILP::modelSCF()
{
	/*
	 * additional (continuous) variables f(i,j,k) represent the amount of "flow" on arc (j;k) by the tour i
	 */

	NumVar3Matrix var_f(env,instance.m);
	for(int i=0; i < instance.m; i++)
	{
		var_f[i] = NumVarMatrix(env, instance.n);
		for(int j=0; j < instance.n; j++)
		{
			var_f[i][j] = IloNumVarArray(env, instance.n);
			for(int k=0; k < instance.n; k++)
			{
				var_f[i][j][k] = IloNumVar(env, Tools::indicesToString( "f_", i, j, k).c_str());
			}
		}
	}

	/*
	 * t(i,j,k) is 1 if the arc from (j,k) is used by the tour i
	 */

	BoolVar3Matrix var_t(env,instance.m);
	IloBoolVarArray var_r(env,instance.m);
	initDecisionVars(var_t,var_r);
	initObjectiveFunction(var_t);
	initConstraints(var_t,var_r);

	/*
	 * Sending out n-1 commodities for every route. n is the number of nodes in this route
	 */

	 for(int i=0;i<instance.m;i++)
	 {
	 	IloExpr myExpr8(env);
	 	IloExpr edgeinExpr(env);
	 	for(int j=1; j < instance.n; j++)
	 	{
	 		myExpr8 += var_f[i][0][j];
	 		for(int k=0; k < instance.n; k++)
	 		{
	 			if(j!=k){
	 				edgeinExpr += var_t[i][k][j] + var_t[i][j][k];
	 			}
	 		}
	 	}
	 	model.add(myExpr8 == (edgeinExpr)/2);
	 	myExpr8.end();
	 	edgeinExpr.end();
	 }

	/*
	 * Leaving one commodity on each node.
	 */

	for(int i=0;i<instance.m;i++)
	{
		for(int j=1;j<instance.n;j++)
		{
			IloExpr incomingExpr(env);
			IloExpr outgoingExpr(env);
			IloExpr edgeinExpr(env);
			for(int k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					incomingExpr += var_f[i][k][j];
					outgoingExpr += var_f[i][j][k];
					edgeinExpr += var_t[i][k][j] + var_t[i][j][k];
				}
			}
			model.add(incomingExpr - outgoingExpr == (edgeinExpr)/2);
			incomingExpr.end();
			outgoingExpr.end();
			edgeinExpr.end();
		}
	}

	/*
	 * the flow must be greater or equal than 0 for all routes
	 */

	for(int i=0;i<instance.m;i++)
	{
		for(int j=0;j<instance.n;j++)
		{
			for(int k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					IloExpr flowExpr(env);
					flowExpr += var_f[i][j][k];
					model.add(flowExpr >= 0);
					flowExpr.end();
				}
			}
		}
	}

	/*
	 * the flow must be smaller than the number of hops
	 */

	for(int i=0;i<instance.m;i++)
	{
		for(int j=0;j<instance.n;j++)
		{
			for(int k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					model.add(var_f[i][j][k] <= instance.n * var_t[i][j][k]);
				}
			}
		}
	}
}

void tcbvrp_ILP::modelMTZ()
{
	/*
	 * additional variables uij , are used to indicate the order in which the nodes are visited on route i
	 */

	NumVarMatrix var_u(env,instance.m);
	for(int i=0; i < instance.m; i++)
	{
		var_u[i] = IloNumVarArray(env, instance.n);
		for(int k=0; k < instance.n; k++)
		{
			var_u[i][k] = IloNumVar(env, Tools::indicesToString( "u_", i, k).c_str());
		}
	}

	/*
	 * t(i,j,k) is 1 if the arc from (j,k) is used by the tour i
	 */

	BoolVar3Matrix var_t(env,instance.m);
	IloBoolVarArray var_r(env,instance.m);
	initDecisionVars(var_t,var_r);
	initObjectiveFunction(var_t);
	initConstraints(var_t,var_r);

	/*
	 * the order must be bigger than 1 and smaller than the number of nodes
	 */

	for(int i=0;i<instance.m;i++)
	{
		IloExpr NumHopsOnRouteIExpr(env);
		for(int j=0;j<instance.n;j++)
		{
			for(int k=0; k < instance.n; k++)
			{
				NumHopsOnRouteIExpr += var_t[i][j][k];
			}
		}
		for(int j=1;j<instance.n;j++)
		{
			IloExpr uSizeExpr(env);
			uSizeExpr += var_u[i][j];
			model.add(uSizeExpr >= 0);
			model.add(uSizeExpr <= NumHopsOnRouteIExpr);
			uSizeExpr.end();
		}
		NumHopsOnRouteIExpr.end();
	}

	/*
	 * if the arc (j,k) on route i is chosen then the order of the node k must be higher than from the node j
	 */

	 for(int i=0;i<instance.m;i++)
	 {
	 	for(int j=1;j<instance.n;j++)
		{
			for(int k=1;(k<instance.n);k++)
			{
				IloExpr uLeftExpr(env);
				IloExpr uRightExpr(env);
				uLeftExpr += var_u[i][j] - var_u[i][k] + 1;
				uRightExpr += (instance.n-1) * (1 - var_t[i][j][k]);
				model.add(uLeftExpr <= uRightExpr);
				uLeftExpr.end();
				uRightExpr.end();
			}
		}
	 }
}

void tcbvrp_ILP::modelMCF()
{
	/*
	 * additional (continuous) variables f(i,j,k) represent the amount of "flow" on arc (j;k) by the tour i
	 */

	 NumVar3Matrix var_f(env,instance.n);
	 for(int i=0; i < instance.n; i++)
	 {
	 	var_f[i] = NumVarMatrix(env, instance.n);
	 	for(int j=0; j < instance.n; j++)
	 	{
	 		var_f[i][j] = IloNumVarArray(env, instance.n);
	 		for(int k=0; k < instance.n; k++)
	 		{
	 			var_f[i][j][k] = IloNumVar(env, Tools::indicesToString( "f_", i, j, k).c_str());
	 		}
	 	}
	 }

	/*
	 * t(i,j,k) is 1 if the arc from (j,k) is used by the tour i
	 */

	BoolVar3Matrix var_t(env,instance.m);
	IloBoolVarArray var_r(env,instance.m);
	initDecisionVars(var_t,var_r);
	initObjectiveFunction(var_t);
	initConstraints(var_t,var_r);

	/*
	 * Sending out 1 commoditie for every used node
	 */

	 for(int k=1; k < instance.n; k++)
	 {
	 	IloExpr myflowExpr(env);
	 	IloExpr edgeinExpr(env);
	 	for(int i=0;i<instance.n;i++)
	 	{
	 		if(i!=k)
	 		{
	 			myflowExpr += var_f[k][i][k];
	 		}
	 		for(int j=0;j<instance.m;j++)
	 		{
	 			edgeinExpr += var_t[j][i][k];
	 		}
	 	}
	 	model.add(myflowExpr == edgeinExpr);
	 	myflowExpr.end();
	 	edgeinExpr.end();
	 }

	/*
	 * assign one comodity to every used node
	 */

	 for(int k=1; k < instance.n; k++)
	 {
	 	IloExpr incomingExpr(env);
	 	IloExpr outgoingExpr(env);
	 	IloExpr edgeinExpr(env);
	 	for(int j=1;j<instance.n;j++)
	 	{
	 		incomingExpr += var_f[k][0][j];
	 		outgoingExpr += var_f[k][j][0];
	 		for(int i=0;i<instance.m;i++)
	 		{
	 			edgeinExpr += var_t[i][j][k];
	 		}
	 	}
	 	model.add(incomingExpr - outgoingExpr == edgeinExpr);
	 	incomingExpr.end();
	 	outgoingExpr.end();
	 	edgeinExpr.end();
	 }

	/*
	 * every node takes the commodody assigned to it.
	 */

	 for(int k=1; k < instance.n; k++)
	 {
	 	for(int j=1;j<instance.n;j++)
	 	{
	 		if(j!=k)
	 		{
	 			IloExpr incomingExpr(env);
	 			IloExpr outgoingExpr(env);
	 			for(int i=0;i<instance.n;i++)
	 			{
	 				if(i!=j)
	 				{
	 					incomingExpr += var_f[k][i][j];
	 					outgoingExpr += var_f[k][j][i];
	 				}
	 			}
	 			model.add(incomingExpr - outgoingExpr == 0);
	 			incomingExpr.end();
	 			outgoingExpr.end();
	 		}
	 	}
	 }

	/*
	 * the flow must be greater or equal than 0 for all routes
	 */

	 for(int k=1; k < instance.n; k++)
	 {
	 	for(int i=0;i<instance.n;i++)
	 	{
	 		for(int j=0;j<instance.n;j++)
	 		{
	 			if(j!=i)
	 			{
	 				IloExpr flowExpr(env);
	 				flowExpr += var_f[k][i][j];
	 				model.add(flowExpr >= 0);
	 				flowExpr.end();
	 			}
	 		}
	 	}
	 }

	/*
	 * the flow must be zero if the node is not used and 1 otherwise
	 */

	 for(int k=1; k < instance.n; k++)
	 {
	 	for(int i=0;i<instance.n;i++)
	 	{
	 		for(int j=0;j<instance.n;j++)
	 		{
	 			if(j!=i)
	 			{
	 				IloExpr flowExpr(env);
	 				IloExpr conExpr(env);
	 				flowExpr += var_f[k][i][j];
	 				for(int l=0;l<instance.m;l++)
	 				{
	 					conExpr += var_t[l][i][j];
	 				}
	 				model.add(flowExpr <= conExpr);
	 				conExpr.end();
	 				flowExpr.end();
	 			}
	 		}
	 	}
	 }
}
