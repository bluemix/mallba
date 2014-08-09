/*****************************************************************
 *
 * First version of PSO for Mallva v0.0
 * Jos� Manuel Garc�a Nieto 01/12/2006
 * File : Pso.cpp BINARY PSO
 * Definition of the PSO skeleton for Mallba
 * For the problem instance: 
 *
 *****************************************************************
 **  Proyecto Fin de Carrera 2005/2006 ****************************/

#include "PSO.hh"

#include "PSO.pro.cc"
#include "PSO.req.cc"

#include <string>
#include <sstream>
#include <iostream>

using namespace std;

int main (int argc, char** argv)
{
	using skeleton PSO;

	char path[MAX_BUFFER];
	int len;
	int longitud;

	//system("clear");

	get_path(argv[0],path);
	len = strlen(path);
	longitud = MAX_BUFFER - len;

	strcat(path,"Config.cfg");
	ifstream f(path);
	if(!f) show_message(2);

	f.getline(&(path[len]),longitud,'\n');
	ifstream funo(path);
	if(!funo)	show_message(2);

	f.getline(&(path[len]),longitud,'\n');
	ifstream fdos(path);
	if(!fdos) show_message(9);
		
	SetUpParams cfg;
	funo >> cfg;
			
	Problem pbm;
	fdos >> pbm;

	Migration mig(pbm.direction(),cfg.migration_rate());

	Solver_Lan solver(pbm,cfg,mig,argc,argv);
	solver.run();

	if (solver.pid()==0)
	{
		solver.show_state();
		cout << "Solution: " << solver.global_best_solution() << endl << endl;
		cout << "Fitness: " << solver.global_best_solution().fitness();

		f.getline(&(path[len]),longitud,'\n');
	  	ofstream fexit(path);
	  	if(!fexit) show_message(2);
	  	fexit << solver.userstatistics();

		cout << endl << endl << " :( ---------------------- THE END --------------- :) " << endl;
	}
	return(0);
}
/*# MainLan */
