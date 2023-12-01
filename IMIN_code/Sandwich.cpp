#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument {
public:
    unsigned int k;
    string dataset;
    string res;
    string algo;
    int Rumor_num;
    double delta;
    double gamma;
    double beta;
    double epsilon;
};

#include "graph.h"
#include "infgraph.h"
#include "Sandwich.h"

void OutputSeedSetToFile(vector<int> seed_set, const Argument& arg)
{
    string seedfile = "results/res_" + arg.dataset;
    ofstream of(seedfile, ios::app);
    //of.open(seedfile);
    for (int seed : seed_set)
    {
        of << seed << endl;
    }
    of << endl;
    of.close();
}

void run_with_parameter(InfGraph& g,  Argument& arg)
{
    CP::CP_based(g, arg);

    //OutputSeedSetToFile(g.seedSet, arg);
}
void Run(int argn, char** argv)
{
    Argument arg;

    arg.k = 0;

    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-k"))
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-dataset"))
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-algo"))
            arg.algo = argv[i + 1];
        if (argv[i] == string("-rumorNum"))
            arg.Rumor_num = atoi(argv[i + 1]);
        if (argv[i] == string("-delta"))
            arg.delta = stod(argv[i + 1]);
        if (argv[i] == string("-gamma"))
            arg.gamma = stod(argv[i + 1]);
        if (argv[i] == string("-beta"))
            arg.beta = stod(argv[i + 1]);
        if (argv[i] == string("-epsilon"))
            arg.epsilon = stod(argv[i + 1]);
    }
    ASSERT(arg.dataset != "");
    string temp_name = arg.dataset.substr(arg.dataset.find_last_of("/") + 1);
    arg.res = "results/res_" + temp_name + "_|S|=" + to_string(arg.Rumor_num) + "_K=" + to_string(arg.k) + "_epsilon=" + to_string(arg.epsilon)
        + "_gamma=" + to_string(arg.gamma) + "_beta=" + to_string(arg.beta) + "_algo=" + arg.algo;
    arg.dataset = arg.dataset + "/";
    string graph_file;
    graph_file = arg.dataset + "graph_ic.inf";
    InfGraph g(arg.dataset, graph_file);
    run_with_parameter(g, arg);
}


int main(int argn, char** argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);

    Run(argn, argv);
}

