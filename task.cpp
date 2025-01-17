#include "task.hpp"
#include "util.hpp"

#include <iostream>
#include <limits>
#include <cmath>
#include <utility>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

using namespace std;

void tasks_sorting(std::vector<Task> &tks)
{
	std::cout << "tasks sorting\n";
	int n = tks.size();
	for (int c = 0; c < n - 1; c++)
	{
		cout << "c: " << c << std::endl;
		for (int d = 0; d < n - c - 1; d++)
		{
			cout << "-->d: " << d << std::endl;
			if (tks[d].get_dline() > tks[d + 1].get_dline())
			{
				cout << "-->--> before swap" << std::endl;
				Task tmp = tks[d];
				tks[d] = tks[d + 1];
				tks[d + 1] = tmp;
				cout << "-->--> after swap" << std::endl;
			}
			else cout << "-->--> no swap\n";
		}
	}
	std::cout << "what happens here\n";
}

template<>
int HasUniqueId<Task>::counter = 0;

Task::Task() : wcet(0), dline(0), period(0), wcrt(numeric_limits<double>::max()), bcrt(numeric_limits<double>::min()) {}

Task::Task(const double c, const double d, const double p)
	: wcet(c), dline(d), period(p), wcrt(numeric_limits<double>::max()), bcrt(numeric_limits<double>::min()) {}

Task::Task(const double c_lb, const double c_ub, const double d, const double p)
	: bcet(c_lb), wcet(c_ub), dline(d), period(p), wcrt(numeric_limits<double>::max()), bcrt(numeric_limits<double>::min()) {}


Task::Task(const Task& t)
	: HasUniqueId<Task>(t), wcet(t.wcet), dline(t.dline), period(t.period), wcrt(t.wcrt), bcrt(t.bcrt), bcet(t.bcet) {}

void Task::print() const {
	cout << "task id: " << get_id() << ", (" << bcet << ", " << wcet << ", " << get_dline() << ", " << period << ") ";
	cout << "wcrt = " << wcrt << ", bcrt = " << bcrt;
	cout << endl;
}


void compute_wcrt(Task& tk, const vector<Task>& hps) {

	double c = tk.get_wcet(), p = tk.get_period();

	double rt = c;
	double wcrt = c;
	double ajrt = c;
	double ajrt0 = c;
	double TOL = 0.2;
	int k = 1;

	bool finished = false;

	do {
		do {
			ajrt0 = ajrt;
			double I = 0;
			for (int i = 0; i < hps.size(); i++) {
				I += ceil(ajrt0 / hps[i].get_period())*hps[i].get_wcet();
			}
			ajrt = I + k*c;

		} while (ajrt > ajrt0);

		rt = ajrt - (k - 1)*p;
		
		// WCRT as worst response time for each iteration
		if (rt > wcrt)
			wcrt = rt;
		else
			finished = true;
		
		k++;

	} while (!finished);

	tk.set_wcrt(wcrt);
}

void compute_bcrt(Task& tk, const vector<Task>& hps) {

	double c = tk.bcet;
	double X = tk.get_period();
	while (true)
	{
		double I = c;
		for (int i = 0; i < hps.size(); i++)
		{
			double p = hps[i].get_period();
			I += ceil((X - p) / p) * hps[i].bcet;
		}
		if (I < X) X = I;
		else break;
	}
	tk.set_bcrt(X);
}

void compute_bp(Task& tk, const vector<Task>& hps) {
	double c = tk.get_wcet(), p = tk.get_period();

	double bp = fmax(c, tk.get_wcrt());

	while (true) {
		//std::cout << "bp on growing: " << bp << std::endl;

		double X = 0;

		for (int i = 0; i < hps.size(); i++)
			X += ceil(bp*1.0 / hps[i].get_period()) * hps[i].get_wcet();
		X += ceil(bp*1.0 / p) * c;

		if (X == bp) {
			break;
			////double Y = X + 0.1;
			//double Y = X + 0.001;

			//X = ceil(Y*1.0/p) * c;
			//for ( int i = 0; i < hps.size(); i++)
			//  X += ceil(Y*1.0/hps[i].get_period()) * hps[i].get_wcet();

			//if ( X >= Y) bp = X;
			//else break;
		}
		else bp = X;
	}
	tk.set_bp(bp);
}

void compute_wcrt_firstjob(Task& tk, const vector<Task>& hps) {

	double c = tk.get_wcet(), p = tk.get_period();
	double rt = c;
	double wcrt = c;
	double ajrt0 = c;

	do {
		ajrt0 = wcrt;
		double I = 0;
		for (int i = 0; i < hps.size(); i++) {
			I += ceil(ajrt0 / hps[i].get_period())*hps[i].get_wcet();
		}
		wcrt = I + c;

	} while (wcrt > ajrt0);

	tk.set_wcrt(wcrt);
}

double computeIdle(const std::vector<Task> & hps, const int L) {

	unsigned int c = 1;

	do {
		Task tk(c, L, L);
		compute_wcrt_firstjob(tk, hps);

		if (tk.wcrt > L)
			break;
		c++;

	} while (c <= L);

	return c - 1;
}

std::vector<Task> read_a_taskset(const string &fname, const int n)
{
	ifstream input(fname);
	vector<Task> tasks;
	string line;

	for (int i = 0; i < n; i++) {
		getline(input, line); // read a line from the input file
		stringstream linestream(line);
		string data;
		getline(linestream, data, '\n');
		vector<string> tokens = split(data, ' ');

		double c, d, p, j;
		c = atof(tokens[0].c_str());
		d = atof(tokens[1].c_str());
		p = atof(tokens[2].c_str());

		if (tokens.size() > 3) // 4th element is jitter
			j = atof(tokens[3].c_str());
		else
			j = 0;

		Task t(c, c, d, p);
		t.set_jitter(j);
		tasks.push_back(t);
	}

	return tasks;
}
