//#include <unistd.h>
//#include <io.h>
#include "task.hpp"
#include "wa.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>

#include <iostream>
using namespace std;



int main (int argc,char *argv[])
{
  // Weakly-hard parameters: to be modified for different tests
  bool job_kill=false;
  unsigned int m = 1;
  unsigned int K = 4;
  
  // Task set parameters
  int  n = 4;
  string fname="example/taskset.dat";

  vector<Task> tasks=read_a_taskset(fname, n);

  for ( int x = 1; x <= n; x++) {
    vector<Task> hps;
    for ( int y = 1; y < x; y++) 
      hps.push_back( tasks[y-1]);
    compute_wcrt(tasks[x-1], hps);
    compute_bp(tasks[x-1], hps);
    compute_bcrt(tasks[x-1], hps);

    tasks[x-1].print();
    if(tasks[x-1].wcrt<=tasks[x-1].dline) continue;
	
	// else: start weakly-hard analysis
    std::cout << "+++++++++++++++++++++++++++++++\n";
    std::cout << "Task " << tasks[x-1].get_id() << " is un-schedulable\n";

    statWA swa;

	std::cout << "Running weakly-hard analysis: " << ", K=" << K << " ... " << std::endl;
	std::cout << "job_kill: " << job_kill << endl;
	
	if (job_kill) {
		swa = wanalysis_kill(tasks[x - 1], hps, K);
		// swa = wanalysis_kill_fast(tasks[x - 1], hps, K); // version without constraints of refined analysis
	}
	else { // job_continue
		swa = wanalysis(tasks[x - 1], hps, m, K);
		// swa = wanalysis_fast(tasks[x - 1], hps, m, K); // version without constraints of refined analysis
	}
	
	std::cout << "*************** Weakly Hard analysis completed with parameters " << swa.m << "," << K << "******\n";
	std::cout << "Time = " << swa.time << " sec " << std::endl;
	std::cout << "+++++++++++++++++++++++++++++++\n";
	std::cout << std::endl << std::endl;
	
  }

  system("pause");
  return 0;
}
