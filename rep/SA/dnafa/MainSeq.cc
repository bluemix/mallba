#include "SA.hh"
#include "ContigBuilder.h"
#include <iostream>
#include <fstream>

extern int Global_cutoff;

int main (int argc, char** argv)
{
	using skeleton SA;

	system("clear");

	if(argc < 4)
		show_message(1);

	ifstream f1(argv[1]);
	if (!f1) show_message(11);

	ifstream f2(argv[2]);
	if (!f2) show_message(12);

	Problem pbm;
	f2 >> pbm;

//cout << pbm;
	Global_cutoff = pbm.cutoff;

	SetUpParams cfg;
	f1 >> cfg;

	Solver_Seq solver(pbm,cfg);
	solver.run();

	if (solver.pid()==0)
	{
		solver.show_state();
		cout << solver.global_best_solution()
		     << " Fitness: " << solver.global_best_solution().fitness() << endl;
		cout << "\n\n :( ---------------------- THE END --------------- :) ";

		ofstream fexit(argv[3]);
		if(!fexit) show_message(13);
		fexit << solver.userstatistics();

	}
	return(0);
}
