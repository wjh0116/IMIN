#pragma once
#pragma once
#include <chrono>
#include <ctime>
#include <ratio>
#include <queue>
#include <ctime>
//#include "infgraph.h"
#define e exp(1)
#define c 2*(exp(1)-2)

using namespace std::chrono;

class Math {
public:
	static double log2(int n) {
		return log(n) / log(2);
	}
	static double logcnk(int n, int k) {
		double ans = 0;
		for (int i = n - k + 1; i <= n; i++)
		{
			ans += log(i);
		}
		for (int i = 1; i <= k; i++)
		{
			ans -= log(i);
		}
		return ans;
	}
};

class CP
{
private:


public:
	

	static void CP_based(InfGraph& g, Argument& arg)
	{

		sfmt_t sfmtSeed;
		sfmt_init_gen_rand(&sfmtSeed, rand());

		double total_spread = 0;
		double total_time = 0;
		string rumor_file;
		rumor_file = arg.dataset + "rumorSet_" + to_string(arg.Rumor_num) + ".txt";
		ifstream inFile(rumor_file);
		if (inFile.is_open())
		{
			int number;
			while (inFile >> number)
			{
				g.rumorSet.push_back(number);
				g.isRumor[number] = true;
			}
			g.isRumor[g.n] = true;
			inFile.close();
		}
		else
		{
			std::cerr << "Unable to open input file!" << std::endl;
		}
		vector<int>CB;//out_neighbor set
		vector<int>is_ON(g.n, 0);
		for (int node : g.rumorSet)
		{
			for (int i = 0; i < g.gT[node].size(); i++)
			{
				int v = g.gT[node][i];
				if (!is_ON[v])
				{
					CB.push_back(v);
					is_ON[v] = 1;
				}
			}
		}
		ofstream of(arg.res, ios::app);
		g.init_hyper_graph();
		double inf1 = g.MC_based_estimate(g.rumorSet, 10000);
		cout << "estimate influence by Mente-Carlo is: " << inf1 << endl;
		g.get_reachable_node(g.Rnode);
		cout << "reachable node size: " << g.Rnode.size() << endl;
		high_resolution_clock::time_point startTime = high_resolution_clock::now();
		int CPnum = 10000;
		if (arg.algo == "SandIMIN")
		{
			if (CB.size() <= arg.k)
			{
				g.seedSet = CB;
			}
			else
			{
				g.seedSet.clear();
				high_resolution_clock::time_point startTime1 = high_resolution_clock::now();
				double inf = g.estimate_inf_byStop(g.rumorSet, g.seedSet, arg.beta, 1.0 / g.n / 6);
				high_resolution_clock::time_point endTime1 = high_resolution_clock::now();
				duration<double> interval = duration_cast<duration<double>>(endTime1 - startTime1);
				cout << "estimate influence time:" << interval.count() << endl;

				cout << "estimate influence by StopAlgorithm is: " << inf << endl;
				double OPT_LB = g.calculate_OPT_lower(arg.k, CB);

				g.opimc_sandwich(arg.k, arg.epsilon, 1.0 / g.n, arg, inf, OPT_LB, 0, 0);
				g.deg_based_heuristic(arg.k, CB);
				g.reset_pro();

				double inf_upper = g.estimate_inf_byStop(g.rumorSet, g.UB_seedSet, arg.gamma, 1.0 / g.n);
				cout << "upper's influence:" << inf_upper << endl;

				double inf_lower = g.estimate_inf_byStop(g.rumorSet, g.LB_seedSet, arg.gamma, 1.0 / g.n);
				cout << "lower's influence:" << inf_lower << endl;

				double inf_or = g.estimate_inf_byStop(g.rumorSet, g.Or_seedSet, arg.gamma, 1.0 / g.n);
				cout << "orginal's influence:" << inf_or << endl;
				if (inf_upper <= inf_lower && inf_upper <= inf_or)
					g.seedSet = g.UB_seedSet;
				else if (inf_lower <= inf_upper && inf_lower <= inf_or)
					g.seedSet = g.LB_seedSet;
				else if (inf_or <= inf_lower && inf_or <= inf_upper)
					g.seedSet = g.Or_seedSet;
			}
		}
		else if (arg.algo == "SandIMIN-")
		{
			if (CB.size() <= arg.k)
			{
				g.seedSet = CB;
			}
			else
			{
				g.seedSet.clear();
				high_resolution_clock::time_point startTime1 = high_resolution_clock::now();
				double inf = g.estimate_inf_byStop(g.rumorSet, g.seedSet, arg.gamma, 1.0 / g.n / 6);
				high_resolution_clock::time_point endTime1 = high_resolution_clock::now();
				duration<double> interval = duration_cast<duration<double>>(endTime1 - startTime1);
				cout << "estimate influence time:" << interval.count() << endl;
				cout << "estimate influence by StopAlgorithm is: " << inf << endl;
				double OPT_LB = g.calculate_OPT_lower(arg.k, CB);

				g.opimc_sandwich(arg.k, arg.epsilon, 1.0 / g.n, arg, inf, OPT_LB, 1, 0);
				
				g.deg_based_heuristic(arg.k, CB);
				g.reset_pro();

				double inf_lower = g.estimate_inf_byStop(g.rumorSet, g.LB_seedSet, arg.beta, 1.0 / g.n);
				cout << "lower's influence:" << inf_lower << endl;

				double inf_or = g.estimate_inf_byStop(g.rumorSet, g.Or_seedSet, arg.beta, 1.0 / g.n);
				cout << "orginal's influence:" << inf_or << endl;
				if (inf_lower <= inf_or)
					g.seedSet = g.LB_seedSet;
				else
					g.seedSet = g.Or_seedSet;
			}
		}
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
		total_time += (double)interval.count();
		cout << "time:" << interval.count() << endl;

		g.reset_pro();
		g.Delete_Node(g.seedSet);
		double inf2 = g.MC_based_estimate(g.rumorSet, 100000);
		double inf3 = inf1 - inf2;
		of << inf1 << "\t" <<  inf2 << "\t" << inf3 << "\t" << (double)interval.count() <<  endl;

		of.close();
	}
};
