#ifndef __WA_HPP__
#define __WA_HPP__
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "task.hpp"
#include <vector>

struct statWA {
  double time = 0;
  bool satisfies_mk = true;
  bool bounded = true;
  int m = 0;
};


std::string  convert_to_string(const int Number);

statWA wanalysis(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

statWA wanalysis_fast(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

statWA wanalysis_kill(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

statWA wanalysis_kill_fast(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

//void wsa(std::vector<Task> &tasks, m, K);

#endif
