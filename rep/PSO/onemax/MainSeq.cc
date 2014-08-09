#include "PSO.hh"

int main (int argc, char** argv)
{
    using skeleton PSO;

    system("clear");

    if(argc < 4)
        show_message(1);

    cout << "argv[1]: " << argv[1] << endl;
    ifstream f1(argv[1]);
    if (!f1) show_message(11);

    cout << "f1: " << f1 << endl;
    /*char buffer[MAX_BUFFER];
    char command[50];

    while(f1.getline(buffer, MAX_BUFFER, '\n'))
    {
        sscanf(buffer," %s ",command);
        cout << "command: " << command << endl;
    }*/


    ifstream f2(argv[2]);
    if (!f2) show_message(12);

    Problem pbm;
    f2 >> pbm;

    ifstream f1_1(argv[1]);
    cout << "f1_1: " << f1_1 << endl;
    SetUpParams cfg;
    f1_1 >> cfg;
    cout << "cfg: " << cfg << endl;

    Solver_Seq solver(pbm,cfg);
    solver.run(); /* number of evaluations */

    if (solver.pid()==0)
    {
        solver.show_state();
        cout << "Solution: " << solver.global_best_solution() << endl
             << " Fitness: " << solver.global_best_solution().fitness() << endl;
        cout << "\n\n :( ---------------------- THE END --------------- :) ";
        cout << endl;
        ofstream fexit(argv[3]);
        if(!fexit) show_message(13);
        fexit << solver.userstatistics();

    }
    return(0);
}
