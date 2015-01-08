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

void tcbvrp_ILP::modelSCF()
{

	// <decision variables>
	
	int i,j,k;	// indices

	typedef IloArray<IloBoolVarArray>	BoolVarMatrix;
	typedef IloArray<BoolVarMatrix>		BoolVar3Matrix;		//3D Bool array

	typedef IloArray<IloNumVarArray>	NumVarMatrix;
	typedef IloArray<NumVarMatrix>		NumVar3Matrix;		//3D Numerical array

	/*
	 * t(i,j,k) is 1 if the arc from (j,k) is used by the tour i
	 */

	BoolVar3Matrix var_t(env,instance.m);
	for(i=0; i < instance.m; i++)
	{
		var_t[i] = BoolVarMatrix(env, instance.n);
		for(j=0; j < instance.n; j++)
		{
			var_t[i][j] = IloBoolVarArray(env, instance.n);
			for(k=0; k < instance.n; k++)
			{
				var_t[i][j][k] = IloBoolVar(env, Tools::indicesToString( "t_", i, j, k).c_str());
			}
		}
	}


	/*
	 * additional (continuous) variables f(i,j,k) represent the amount of "flow" on arc (j;k) by the tour i
	 */
	NumVar3Matrix var_f(env,instance.m);
	for(i=0; i < instance.m; i++)
	{
		var_f[i] = NumVarMatrix(env, instance.n);
		for(j=0; j < instance.n; j++)
		{
			var_f[i][j] = IloNumVarArray(env, instance.n);
			for(k=0; k < instance.n; k++)
			{
				var_f[i][j][k] = IloNumVar(env, Tools::indicesToString( "f_", i, j, k).c_str());
			}
		}
	}

	// </decision variables>

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

	// <objective function>
	IloExpr objFunction(env);
	for(i=0;i<instance.m;i++)
	{
		for(j=0; j< instance.n; j++)
		{
			for(k=0; k< instance.n; k++)
			{
				objFunction += var_t[i][j][k] * instance.getDistance(j, k);
			}
		}
	}

	model.add(IloMinimize(env, objFunction));
	objFunction.end();
	// </objective function>

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

	// <constraints>

	/*
	 * Each demand node has to have an outgoing arc which goes to a supply node or the originator
	 */
	
	for(j=1; j< instance.n; j++)
	{
		if(instance.isDemandNode(j))
		{
			IloExpr myExpr(env);
			for(i=0;i<instance.m;i++)
			{
				for(k=0; k < instance.n; k++)
				{
					if(instance.isSupplyNode(k) || k == 0)
					{
						myExpr += var_t[i][j][k];
					}
				}
			}
			model.add(myExpr == 1);
			myExpr.end();
		}
	}

	/*
	 * Each node is not allowed to go to itself
	 */
	
	// for(j=0; j< instance.n; j++)
	// {
	// 	if(instance.isDemandNode(j))
	// 	{
	// 		IloExpr noSelfExpr(env);
	// 		for(i=0;i<instance.m;i++)
	// 		{
	// 			for(k=1; k < instance.n; k++)
	// 			{
	// 					noSelfExpr += var_t[i][j][k];
	// 			}
	// 		}
	// 		model.add(noSelfExpr == 0);
	// 		noSelfExpr.end();
	// 	}
	// }

	/*
	 * Each supply node can only go to at most one demand node
	 */

	for(j=1; j< instance.n; j++)
	{
		if(instance.isSupplyNode(j))
		{
			IloExpr myExpr1(env);
			for(i=0;i<instance.m;i++)
			{
				for(k=0; k < instance.n; k++)
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
	 * A supply node is not allowed to go to an other supply node or the originator
	 */

	for(j=1; j< instance.n; j++)
	{
		if(instance.isSupplyNode(j))
		{
			IloExpr myExpr2(env);
			for(i=0;i<instance.m;i++)
			{
				for(k=0; k < instance.n; k++)
				{
					if(instance.isSupplyNode(k) || k==0)
					{
						myExpr2 += var_t[i][j][k];
					}
				}
			}
			model.add(myExpr2 == 0);
			myExpr2.end();
		}
	}

	/*
	 * The originator is not allowed to have more than m outgoing arcs to supply nodes
	 */

	IloExpr myExpr3(env);
	for(i=0;i<instance.m;i++)
	{
		for(k=0; k < instance.n; k++)
		{
			if(instance.isSupplyNode(k))
			{
				myExpr3 += var_t[i][0][k];
			}
		}
	}
	model.add(myExpr3 <= instance.m);
	myExpr3.end();

	/*
	 * If the ingoing arc in one node is from a tour the outgoing arc has to be from the same tour
	 */

	for(j=0; j< instance.n; j++)
	{
		for(i=0;i<instance.m;i++)
		{
			IloExpr myExpr4(env);
			IloExpr myExpr5(env);
			for(k=0; k < instance.n; k++)
			{
				myExpr4 += var_t[i][k][j];

			}

			for(k=0; k < instance.n; k++)
			{
				myExpr5 += var_t[i][j][k];
			}

			model.add(myExpr4 == myExpr5);
			myExpr4.end();
			myExpr5.end();
		}
	}

	/*
	 * The originator is not allowed to go to a demand node or the originator itself
	 */

	IloExpr originNotToDemandExpr(env);
	for(i=0;i<instance.m;i++)
	{
		for(k=0; k < instance.n; k++)
		{
			if(instance.isDemandNode(k) || k == 0)
			{
				originNotToDemandExpr += var_t[i][0][k];
			}

		}

	}
	model.add(originNotToDemandExpr == 0);
	originNotToDemandExpr.end();

	/*
	 * A tour must be finished under the maximum time
	 */

	for(i=0;i<instance.m;i++)
	{
		IloExpr maxTimeExpr(env);
		for(j=0;j<instance.n;j++)
		{
			for(k=0; k < instance.n; k++)
			{
				maxTimeExpr += var_t[i][j][k] * instance.getDistance(j, k);
			}

		}
		model.add(maxTimeExpr <= instance.T);
		maxTimeExpr.end();
	}

	//<SCF>

	/*
	 * Sending out n-1 commodities for every route. n is the number of nodes in this route
	 */

	for(i=0;i<instance.m;i++)
	{
		IloExpr myExpr8(env);
		IloExpr myExpr9(env);
		for(j=1;j<instance.n;j++)
		{
			myExpr8 += var_f[i][0][j];
		}
		for(j=0;j<instance.n;j++)
		{
			for(k=0; k < instance.n; k++)
			{
				myExpr9 += var_t[i][j][k];
			}

		}
		model.add(myExpr8 == myExpr9 - 1);
		myExpr8.end();
		myExpr9.end();
	}

	/*
	 * Leaving one commodity on each node.
	 */

	for(i=0;i<instance.m;i++)
	{
		for(j=1;j<instance.n;j++)
		{
			IloExpr myExpr10(env);
			IloExpr myExpr11(env);
			for(k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					myExpr10 += var_f[i][k][j];
				}

			}
			for(k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					myExpr11 += var_f[i][j][k];
				}

			}
			model.add(myExpr10 - myExpr11 == 1);
			myExpr10.end();
			myExpr11.end();


		}
	}

	/*
	 * the flow must be greater than 0 for all routes
	 */

	for(i=0;i<instance.m;i++)
	{
		for(j=0;j<instance.n;j++)
		{
			for(k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					IloExpr myExpr12(env);
					myExpr12 += var_f[i][j][k];

					model.add(myExpr12 >= 0);
					myExpr12.end();
				}

			}

		}
	}

	/*
	 * the flow must be smaller than the number of hops in this route (ERROR!)
	 */

	for(i=0;i<instance.m;i++)
	{
		IloExpr myExpr13(env);

		for(j=0;j<instance.n;j++)
		{
			for(k=0; k < instance.n; k++)
			{
				myExpr13 += var_t[i][j][k];

			}

		}

		for(j=0;j<instance.n;j++)
		{
			for(k=0; k < instance.n; k++)
			{
				if(j!=k)
				{
					IloExpr myExpr14(env);
					myExpr14 += var_f[i][j][k];

					model.add(myExpr14 <= (myExpr13-1));
					myExpr14.end();
				}

			}

		}
	}
	// </SCF>
	// </constraints>
}

void tcbvrp_ILP::modelMTZ()
{
	// ++++++++++++++++++++++++++++++++++++++++++
	// TODO build Miller-Tucker-Zemlin model
	// ++++++++++++++++++++++++++++++++++++++++++
}

void tcbvrp_ILP::modelMCF()
{
	// ++++++++++++++++++++++++++++++++++++++++++
	// TODO build multi commodity flow model
	// ++++++++++++++++++++++++++++++++++++++++++
}



