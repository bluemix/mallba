/*****************************************************************
 *
 * First version of PSO for Mallva v0.0
 * José Manuel García Nieto 01/12/2006
 * File : Pso.cpp BINARY PSO
 * Definition of the PSO skeleton for Mallba
 * For the problem instance:
 *
 *****************************************************************
 **  Proyecto Fin de Carrera 2005/2006 ****************************/

#include "PSO.hh"

int main (int argc, char** argv)
{
    using skeleton PSO;

    char path[MAX_BUFFER];
    int len;
    int longitud;

    //system("clear");

    //char the_path[256];

    getcwd(path, MAX_BUFFER);
    strcat(path, "/");
    //strcat(path, argv[0]);

    //printf("%s\n", the_path);

    //cout << "argv[0]: " << argv[0] << endl;

    //get_path(argv[0],path);

    len = strlen(path);
    longitud = MAX_BUFFER - len;

    cout << "longitud: " << longitud << endl;

    strcat(path,"Config.cfg");
    //system("pwd");
    cout << "path: " << path << endl;

    ifstream f(path);
    if(!f) show_message(2);

    f.getline(&(path[len]),longitud,'\n');
    ifstream f1(path);
    if(!f1)	show_message(2);

    f.getline(&(path[len]),longitud,'\n');
    
    ifstream f2(path);
    if(!f2) show_message(9);

    SetUpParams cfg;
    f1 >> cfg;
    //cout << cfg << endl;
    Problem pbm;
    f2 >> pbm;


    Migration mig(pbm.direction(),cfg.migration_rate());

    Solver_Lan solver(pbm,cfg,mig,argc,argv);
    solver.run();

    if (solver.pid()==0)
    {
        solver.show_state();
        cout << "Solution: " << solver.global_best_solution() << " Fitness: " << solver.global_best_solution().fitness();

        f.getline(&(path[len]),longitud,'\n');
        ofstream fexit(path);
        if(!fexit) show_message(2);
        fexit << solver.userstatistics();

        cout << endl << endl << " :( ---------------------- THE END --------------- :) " << endl;
    }
    return(0);
}
/*# MainLan */
