#ifndef __TASK_HPP__
#define __TASK_HPP__

#include <string>
#include <vector>

#include "has_unique_id.hpp"


class Task : public HasUniqueId<Task> {

  private:
    std::string name;

  public:

    double wcet, bcet, dline, period;
    double wcrt, bcrt;
    double bp;
    double end_r;
	double jitter = 0;


    Task ();
    Task (const double c, const double d, const double p); 
    Task (const double c_lb, const double c_ub, const double d, const double p);
    Task (const Task& t);


    void set_wcet(const double w) { wcet = w;}
    void set_wcrt(const double w) { wcrt = w;}
    void set_bcrt(const double b) { bcrt = b;}
    void set_bp(const double b) { bp = b;}
    void set_end_r(const double b) { end_r = b;}
	void set_jitter(const double j) { jitter = j;}

    double get_wcet() const { return wcet;}
    double get_bcet() const { return bcet;}
    double get_dline() const { return dline;}
    double get_period() const { return period;}
    double get_wcrt() const { return wcrt;}
    double get_bcrt() const { return bcrt;}
    double get_bp() const { return bp;}
    double get_end_r() const { return end_r;}
	double get_jitter() const { return jitter; }

    void print() const ;

};


void compute_bp(Task& tk, const std::vector<Task>& hps);
void compute_wcrt(Task& tk, const std::vector<Task>& hps);
void compute_bcrt(Task& tk, const std::vector<Task>& hps);

void tasks_sorting(std::vector<Task> &tks);

std::vector<Task> read_a_taskset(const std::string &fname, const int n);

#endif
