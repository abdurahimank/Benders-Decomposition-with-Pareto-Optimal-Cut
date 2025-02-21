/************************************************************************************
 * This example on Benders Decomposition is taken from the tutorial on Benders Decomposition.
 * Refer to Benders-Slides.pdf
 * Minimize 10Y + 18X1 + 8X2 + 20X3
 * s.t.:
 * 3X1 + X2 + X3 + 2Y >= 12
 * X1 + X2 + 4X3 + Y >= 10
 * X1, X2, X3 >= 0, Y integer in {2,12}
 ***********************************************************************************/
#pragma warning(disable : 4996)//For Visual Studio 2012
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>//for vectors
#include <math.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

using namespace std;

ILOSTLBEGIN

int main(int argc, char** argv)
{
	IloEnv env;
	try
	{
		////////DECISION VARIABLES AND PARAMETERS FOR SUBPROBLEM///////////////
		IloNumVarArray X_dual(env, 2, 0, IloInfinity, ILOFLOAT);
		IloNum Y_val;
		IloNum theta_val;
		IloNumVarArray aux_dual(env, 2, 0, IloInfinity, ILOFLOAT); // new

		////////DECISION VARIABLES AND PARAMETERS FOR MASTER PROBLEM///////////
		IloNumVar Y(env, 2, 12, ILOINT);
		IloNumArray X_dual_val(env, 2);
		//IloNumVar theta_var(env, -IloInfinity, IloInfinity, ILOFLOAT);
		IloNumVar theta_var(env, 0, IloInfinity, ILOFLOAT);
		IloNumArray D_aux_val(env, 2, 0, IloInfinity); // new

		//////////DEVELOP GENERIC MODEL //////////////////////////

		//////SET MASTER PROBLEM///////////////////////////
		IloModel model_master(env);
		IloExpr Objective_master(env);
		Objective_master = theta_var + 10 * Y;
		model_master.add(IloMinimize(env, Objective_master));
		IloCplex cplex_master(env);
		Objective_master.end();
		cplex_master.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_master.setWarning(env.getNullStream()); //This is to supress warning messages on screen

		/////////SET SUBPROBLEM (DUAL FORMULATION)//////////////////////
		IloModel model_sub(env);
		IloObjective Objective_sub = IloMaximize(env);
		model_sub.add(Objective_sub);

		model_sub.add(3 * X_dual[0] + X_dual[1] <= 18);
		model_sub.add(X_dual[0] + X_dual[1] <= 8);
		model_sub.add(X_dual[0] + 4 * X_dual[1] <= 20);


		IloCplex cplex_sub(model_sub);
		IloNum eps = cplex_sub.getParam(IloCplex::EpInt);//Integer tolerance for MIP models; 
		//default value of EpInt remains 1e-5 http://www.iro.umontreal.ca/~gendron/IFT6551/CPLEX/HTML/relnotescplex/relnotescplex12.html
		cplex_sub.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_sub.setWarning(env.getNullStream()); //This is to supress warning messages on screen


		/////////Auxiliary Subproblem (For Pareto-Optimal Cut)//////////////////////
		IloModel model_aux_sub(env);
		IloExpr obj_aux_sub(env);

		// Y should be in the relative interior (2 < Y < 12), then we will have only one solution
		// choosing y as 5 (obj fun of sub problem = max (12-y)u1 + (10-y)u2
		obj_aux_sub = 7 * aux_dual[0] + 5 * aux_dual[1];
		model_aux_sub.add(IloMaximize(env, obj_aux_sub));
		model_aux_sub.add(3 * aux_dual[0] + aux_dual[1] <= 18);
		model_aux_sub.add(aux_dual[0] + aux_dual[1] <= 8);
		model_aux_sub.add(aux_dual[0] + 4 * aux_dual[1] <= 20);

		IloCplex cplex_aux_sub(env);
		obj_aux_sub.end();
		cplex_aux_sub.setOut(env.getNullStream());
		cplex_aux_sub.setWarning(env.getNullStream());



		/////////BEGIN ITERATIONS/////////////////////////////////
		IloNum GAP = IloInfinity;
		theta_val = 0;
		Y_val = 2;
		IloNum sub_obj_val = 0;
		IloNum Upper_bound = IloInfinity;
		IloNum Lower_bound = 0;
		GAP = Upper_bound - Lower_bound;
		cout << "Y = " << Y_val << endl;
		IloInt Iter = 0;

		//while( Iter < MaxCut )
		while (Upper_bound - Lower_bound > eps)
		{
			Iter++;
			cout << "=========================================" << endl;
			cout << "============ITERATION " << Iter << "==============" << endl;
			//Define Object Function for the Dual of the Sub problem
			IloExpr sub_obj(env);
			sub_obj = (12 - 2 * Y_val) * X_dual[0] + (10 - Y_val) * X_dual[1];
			Objective_sub.setExpr(IloMaximize(env, sub_obj));
			//cplex_sub.setParam(cplex_sub.PreInd, 0);   //Disable presolve, otherwise, if dual is infeasible, 
																   //we don't know if prime is unbounded or infeasible
			cout << "SOLVING SUB PROBLEM" << endl;
			cplex_sub.solve();
			cout << "Sub Problem Solution Status: " << cplex_sub.getCplexStatus() << endl;
			if (cplex_sub.getCplexStatus() == CPX_STAT_OPTIMAL)
			{// Dual subproblem is bounded; Add Optimality Cut to the Master Problem
				cplex_sub.getValues(X_dual_val, X_dual);
				cout << "X_dual = " << X_dual_val << endl;
				sub_obj_val = cplex_sub.getObjValue();
				cout << "sub_obj_val = " << sub_obj_val << endl;
				Upper_bound = IloMin(Upper_bound, (10 * Y_val + sub_obj_val));
				cout << "Upper_bound = " << Upper_bound << endl;



				cplex_aux_sub.setParam(IloCplex::PreInd, IloFalse);  // Disable presolve
				cplex_aux_sub.setParam(IloCplex::RootAlg, IloCplex::Barrier);  // Try a different solver


				///// Solve Auxiliary Subproblem for Pareto-Optimal Cut
				cplex_aux_sub.extract(model_aux_sub);
				IloRange r1(env, sub_obj_val, (12 - 2 * Y_val) * aux_dual[0] + (10 - Y_val) * aux_dual[1], sub_obj_val);
				IloExpr aux_obj(env);

				model_aux_sub.add(r1);
				cplex_aux_sub.solve();
				cplex_aux_sub.getValues(D_aux_val, aux_dual);


				cout << "SOLVING AUXILIARY SUBPROBLEM" << endl;
				cplex_aux_sub.extract(model_aux_sub);
				cplex_aux_sub.solve();
				cout << "Auxiliary Subproblem Solution Status: " << cplex_aux_sub.getCplexStatus() << endl;

				IloNum cut_coeff = (2 * D_aux_val[0] + D_aux_val[1]);
				IloNum cut_rhs = (12 * D_aux_val[0] + 10 * D_aux_val[1]);
				cout << "Pareto-Optimal Cut Added: theta + " << cut_coeff << " * Y >= " << cut_rhs << endl;
				model_aux_sub.remove(r1);
				model_master.add(theta_var + cut_coeff * Y >= cut_rhs);




				//Add Cut to the Master Problem
				//cout << "Cut Added to Master Problem: " << "theta + " << (2 * X_dual_val[0] + X_dual_val[1]) << " Y >= " << 6 * X_dual_val[0] + 10 * X_dual_val[1] << endl;
				//model_master.add(theta_var + (2 * X_dual_val[0] + X_dual_val[1]) * Y >=
					//6 * X_dual_val[0] + 10 * X_dual_val[1]);
			}
			cout << "SOLVING MASTER PROBLEM" << endl;
			cout << "Master Problem Solution Status: " << cplex_master.getCplexStatus() << endl;
			cplex_master.extract(model_master);
			if (!cplex_master.solve())
			{
				cout << "Failed" << endl;
				throw(-1);
			}
			Y_val = cplex_master.getValue(Y);
			theta_val = cplex_master.getValue(theta_var);
			cout << "theta_var = " << theta_val << endl;
			cout << "Y = " << Y_val << endl;
			Lower_bound = 10 * Y_val + theta_val;
			cout << "Lower_bound = " << Lower_bound << endl;
		}//while(Upper_bound - Lower_bound > eps)
		model_master.end();
		model_sub.end();
		cplex_master.end();
		cplex_sub.end();
	}//try
	catch (IloException& e)
	{
		env.out() << "ERROR: " << e << endl;
	}
	catch (...)
	{
		env.out() << "Unknown exception" << endl;
	}
	env.end();
	return 0;
}