#include "wa.hpp"

//#include <ilcplex/ilocplex.h>
//ILOSTLBEGIN

#include <iostream>
#include <string>
#include <limits>
#include <algorithm>
using namespace std;

/***
* The weakly hard real-time analysis for periodic task systems.
* Given 'K' successive activations of "ti", can the number of deadline
* misses exceed 'm'? When missing its deadline, the job is killed...
*/
statWA wanalysis_kill(const Task& ti, const std::vector<Task>& hps, const int m, const int K)
{

	cout << "Weakly hard real-time analysis of a task ti" << endl;
	double ci = ti.get_wcet(), di = ti.get_dline(), pi = ti.get_period(), Ri = ti.get_wcrt(), bri = ti.get_bcrt();
	double BP = ti.get_bp(); int Ni = ceil(BP / pi);
	double ui = ci / pi;
	double ji = ti.get_jitter();
	//double ji = 0;

	// the vector of wcet, dline, period and utilization of each higher priority task
	vector<double> c, d, p, u, br, jitt;
	for (unsigned int j = 1; j <= hps.size(); j++) {
		double cj = hps[j - 1].get_wcet(), dj = hps[j - 1].get_dline(), pj = hps[j - 1].get_period();
		double brj = hps[j - 1].get_bcrt();
		double jitj = hps[j - 1].get_jitter();
		//double jitj = 0;
		c.push_back(cj); d.push_back(dj); p.push_back(pj); u.push_back(cj / pj);
		br.push_back(brj);
		jitt.push_back(jitj);
	}

	const double M = 10000000, sigma1 = 1, epsilon = 1;

	vector<Task> copy_hps(hps);
	copy_hps.push_back(ti);
	vector<double> idle_lb; // minimum level-i idle time w.r.c time lengths of pi, 2*pi, 3*pi, ..., K*pi
	for (unsigned int i = 1; i <= K; i++) {
		idle_lb.push_back(computeIdle(copy_hps, i*pi));
	}



	IloEnv env;
	statWA swa;

	try {

		IloModel model(env);

		/**
		* ASSUMPTIONS
		*
		* A1: the job in prior to the 1st job is schedulable (otherwise, we can
		* always advance the job windows without reducing the number of deadline
		* misses)
		*
		* A2: the 1st job must misses its deadline, otherwise, we can always
		* postpone the beginning of the problem windows without reducing the
		* number of deadline misses)
		*
		*/


		/*********************** VARIABLES **********************************/

		// (4.1.1) Activation instant
		IloNumVarArray a(env, 0, 0, IloInfinity, ILOFLOAT); 
		for (unsigned int k = 1; k <= K + 1; k++) {
			string name_a = "a" + convert_to_string(k);
			a.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_a.c_str())); 
			if (k > 1)
				model.add(IloRange(env, pi - ji, a[k - 1] - a[k - 2], pi + ji));
		}

		// (4.1.2) Busy windows 
		IloNumVarArray L(env, 0, 0, IloInfinity, ILOFLOAT); 
		for (unsigned int k = 1; k <= K + 1; k++) {
			string name_bp = "bp" + convert_to_string(k);
			L.add(IloNumVar(env, 0, pi - bri + ji, ILOFLOAT, name_bp.c_str())); 
		}
		// First activation
		model.add(IloRange(env, 0, a[0] - L[0], ji));

		// (4.1.3)a offsets 
		IloNumVarArray alpha(env, 0, 0, IloInfinity, ILOFLOAT); 
		for (unsigned int j = 1; j <= hps.size(); j++) {
			string name_alpha = "alpha" + convert_to_string(j);
			alpha.add(IloNumVar(env, 0, p[j - 1] - br[j - 1] + jitt[j - 1], ILOFLOAT, name_alpha.c_str())); //ILOINT
		}


		//// (4.1.3)b check that at least one offset is 0
		IloNumVarArray b_off(env, 0, 0, 1, ILOBOOL);
		IloExpr sum_boff(env);
		for (unsigned int j = 1; j <= hps.size() + 1; j++) {
			string name_b_off = "b_off" + convert_to_string(j);
			b_off.add(IloNumVar(env, 0, 1, ILOBOOL, name_b_off.c_str()));
			sum_boff += b_off[j - 1];
			if (j <= hps.size())
				model.add(IloRange(env, -IloInfinity, alpha[j - 1] - M*b_off[j - 1], 0));
			else
				model.add(IloRange(env, -IloInfinity, L[0] - M*b_off[j - 1], 0));
		}
		model.add(IloRange(env, 0, sum_boff, hps.size()));

		// (4.1.4) idle times 
		IloNumVarArray idle(env, 0, 0, IloInfinity, ILOFLOAT);
		for (unsigned int k = 1; k <= K; k++) {
			string name_idle = "idle" + convert_to_string(k);
			idle.add(IloNumVar(env, 0, pi, ILOFLOAT, name_idle.c_str()));
		}

		// (5.1.1) finish times
		IloNumVarArray f(env, 0, 0, IloInfinity, ILOFLOAT);
		for (unsigned int k = 1; k <= K; k++) {
			string name_ft = "ft" + convert_to_string(k);
			f.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_ft.c_str()));

			model.add(IloRange(env, 0, f[k - 1] - a[k - 1], di));

			//if (k > 1)
				//model.add(IloRange(env, pi - di + bri, f[k - 1] - f[k - 2], pi - bri + di));
		}


		// (5.1.2) extra finish times 
		IloNumVarArray fe(env, 0, 0, IloInfinity, ILOFLOAT);
		for (unsigned int k = 1; k <= K; k++) {
			string name_ft = "fe" + convert_to_string(k);
			fe.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_ft.c_str()));

			model.add(IloRange(env, 0, fe[k - 1] - a[k - 1], Ri));
			model.add(IloRange(env, 0, fe[k - 1] - f[k - 1], Ri));
		}

		// delta-execution
		IloNumVarArray delta(env, 0, 0, IloInfinity, ILOFLOAT); 
		for (unsigned int k = 1; k <= K; k++) {
			string name_ft = "delta" + convert_to_string(k);
			delta.add(IloNumVar(env, 0.0, ci, ILOFLOAT, name_ft.c_str())); 
		}


		// Schedulability 
		// b == 1 --> the job is not schedulable; b == 0 --> schedulable
		IloNumVarArray b(env, 0, 0, 1, ILOBOOL);
		for (unsigned int k = 1; k <= K; k++) {
			string name = "b" + convert_to_string(k);
			b.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));

			IloExpr dk(env); dk += L[0] + (k - 1)*pi + di;
			// job-kill: b=1 --> dk < fe
			model.add(IloRange(env, -IloInfinity, dk - fe[k - 1] - M*(1.0 - 1.0*b[k - 1]) + sigma1, 0));
			// job-kill: b=0 --> delta = ci
			// b=1 ---> delta < ci
			model.add(IloRange(env, -IloInfinity, delta[k - 1] - ci + 1.0*b[k - 1], 0));
			model.add(IloRange(env, 0, delta[k - 1] - ci + b[k - 1] * M, IloInfinity));
			// job-kill: b=0 --> f = fe
			model.add(IloRange(env, 0, f[k - 1] - fe[k - 1] + (b[k - 1])*M, IloInfinity));
			model.add(IloRange(env, -IloInfinity, f[k - 1] - fe[k - 1], 0));
			// b=0 --> f-a_k >= ri
			model.add(IloRange(env, 0, f[k - 1] + (b[k - 1])*M - a[k - 1] - bri, IloInfinity));
			// b[k]=0 --> L[k+1]<=pi-ri
			if (k != K)
				model.add(IloRange(env, -IloInfinity, L[k] + (b[k - 1])*(di - bri) - (pi - bri), 0));
		}
		// b[0] == 1: J_1 misses its deadline
		model.add(IloRange(env, 1, b[0], 1)); 


	    // (4.1.7) interference from the previous job 
		// job-kill: no interference from prior job!


		// (5.1.3) counting jobs from a higher priority task 
		//// nf_j_k: #jobs of tau_j in [0,f_k)
		//// nL_j_k: #jobs of tau_j in: [0, r_k-L_k) if b_{k-1}=0;
		//// nfe_j_k: #jobs of tau_j in [0, f_k^{e}) job-kill
		std::vector<IloNumVarArray> nf, nL, nfe;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			nf.push_back(IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
			nL.push_back(IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
			nfe.push_back(IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
			for (unsigned int k = 1; k <= K; k++) {
				string name_f =
					"nf_" + convert_to_string(j) + "_" + convert_to_string(k);
				nf[j - 1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_f.c_str()));

				string name_b =
					"nL_" + convert_to_string(j) + "_" + convert_to_string(k);
				nL[j - 1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_b.c_str()));

				string name_e =
					"nfe_" + convert_to_string(j) + "_" + convert_to_string(k);
				nfe[j - 1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_e.c_str()));

			}
			// nL_{K+1}
			string name_b = "nL_" + convert_to_string(j) + "_" + convert_to_string(K + 1);
			nL[j - 1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_b.c_str()));
		}

		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K; k++) {
				//// a coarse upper bound (uf) 
				model.add(IloRange(env, -IloInfinity, nf[j - 1][k - 1] - ceil((pi - bri + (k - 1)*pi + ji + di) / p[j - 1]), 0));

				//// nf 
				model.add(IloRange(env, 0, nf[j - 1][k - 1] * p[j - 1] - (f[k - 1] - alpha[j - 1] - jitt[j - 1]), IloInfinity));
				model.add(IloRange(env, -IloInfinity, nf[j - 1][k - 1] * p[j - 1] - (f[k - 1] - alpha[j - 1] + jitt[j - 1]) - p[j - 1] + sigma1, 0));

				// job-kill
				model.add(IloRange(env, -IloInfinity, nfe[j - 1][k - 1] - ceil((pi - bri + (k - 1)*pi + ji + Ri) / p[j - 1]), 0));

				//// nfe
				model.add(IloRange(env, 0, nfe[j - 1][k - 1] - (fe[k - 1] - alpha[j - 1] - jitt[j - 1]) / p[j - 1], IloInfinity));
				model.add(IloRange(env, -IloInfinity, nfe[j - 1][k - 1] * p[j - 1] + sigma1 - (fe[k - 1] - alpha[j - 1] + jitt[j - 1]) - p[j - 1], 0));

				if (k > 1) {
					//// a coarse upper bound (uf) 
					model.add(IloRange(env, -IloInfinity, nL[j - 1][k - 1] - ceil((pi - bri + (k - 1)*pi) / p[j - 1]), 0));

					//// nL [11]
					model.add(IloRange(env, 0, nL[j - 1][k - 1] - (a[k - 1] - L[k - 1] - alpha[j - 1] - jitt[j - 1]) / p[j - 1], IloInfinity));
					model.add(IloRange(env, -IloInfinity, nL[j - 1][k - 1] * p[j - 1] + sigma1 - (a[k - 1] - L[k - 1] - alpha[j - 1] + jitt[j - 1]) - p[j - 1], 0));
				}
				else {
					model.add(IloRange(env, 0, nL[j - 1][k - 1], 0));
				}
			}

			//// a coarse upper bound (uf) 
			model.add(IloRange(env, -IloInfinity, nL[j - 1][K] - ceil((pi - bri + K*pi + ji) / p[j - 1]), 0));

			//// nL 
			model.add(IloRange(env, 0, nL[j - 1][K] - (a[K] - L[K] - alpha[j - 1] - jitt[j - 1]) / p[j - 1], IloInfinity));
			model.add(IloRange(env, -IloInfinity, nL[j - 1][K] * p[j - 1] + sigma1 - (a[K] - L[K] - alpha[j - 1] + jitt[j - 1]) - p[j - 1], 0));
		}


		// (4.1.9) refining the job counting for a higher priority task

		//// nfb_j_k[]
		vector< std::vector<IloNumVarArray> > nfb;
		double ubf = di + pi;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			vector<IloNumVarArray> v;
			for (unsigned int k = 1; k <= K; k++) {
				IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
				for (unsigned int i = 1; i <= ceil(ubf / p[j - 1]); i++) {
					string name = "nfb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
					vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
				}
				v.push_back(vv);
			}
			nfb.push_back(v);
		}
		////// nL_j_k + Delta = nf_j_k
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K; k++) {
				IloExpr Delta(env);
				for (unsigned int i = 1; i <= ceil(ubf / p[j - 1]); i++)
					Delta += nfb[j - 1][k - 1][i - 1];
				model.add(IloRange(env, 0, nL[j - 1][k - 1] + Delta - nf[j - 1][k - 1], 0));
			}
		}
		////// precedent constraint on ""nfb"s
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K; k++) {
				for (unsigned int i = 1; i <= ceil(ubf / p[j - 1]); i++) {
					if (i > 1)
						model.add(IloRange(env, 0, nfb[j - 1][k - 1][i - 1] - nfb[j - 1][k - 1][i - 2], 1));
				}
			}
		}

		//// nLb_j_k[]
		vector< std::vector<IloNumVarArray> > nLb;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			vector<IloNumVarArray> v;
			for (unsigned int k = 1; k <= K + 1; k++) {
				IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
				for (unsigned int i = 1; i <= ceil((pi) / p[j - 1]); i++) {
					string name = "nLb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
					vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
				}
				v.push_back(vv);
			}
			nLb.push_back(v);
		}

		////// nf[k-2]+Delta_p == nL[k-1] 
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 2; k <= K + 1; k++) {
				IloExpr Delta_p(env);
				for (unsigned int i = 1; i <= ceil((pi) / p[j - 1]); i++)
					Delta_p += nLb[j - 1][k - 1][i - 1];
				model.add(IloRange(env, 0, nf[j - 1][k - 2] + Delta_p - nL[j - 1][k - 1], 0));
			}
		}

		////// nLb[j][k][i] >= nLb[j][k][i-1]
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K + 1; k++) {
				for (unsigned int i = 1; i <= ceil((pi) / p[j - 1]); i++) {
					if (i > 1)
						model.add(IloRange(env, 0, nLb[j - 1][k - 1][i - 1] - nLb[j - 1][k - 1][i - 2], 1));
				}
			}
		}

		// ***************** REFINED ANALYSIS ****************************** //
		// (5.2.1): "sjk" is the activation time of the last job of task j before fe_k

		std::vector<IloNumVarArray> sjk;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			sjk.push_back(IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
			for (unsigned int k = 1; k <= K; k++) {
				string name_s = "sjk_" + convert_to_string(j) + "_" + convert_to_string(k);
				sjk[j - 1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_s.c_str()));
			}
		}

		// n_s_j_k is the number of activations of task j in the interval [0, sjk)
		std::vector<std::vector<IloNumVarArray>>nsjk;
		nsjk.resize(hps.size());
		for (unsigned int s = 1; s <= hps.size(); s++) {
			for (unsigned int j = 1; j <= hps.size(); j++) {
				nsjk[s - 1].push_back(IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
				for (unsigned int k = 1; k <= K; k++) {
					string name = "nsjk_" + convert_to_string(s) + convert_to_string(j) + convert_to_string(k);
					nsjk[s - 1][j - 1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name.c_str()));
				}
			}
		}

		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K; k++) {
				// sjk definition
				model.add(IloRange(env, 0, sjk[j - 1][k - 1] - alpha[j - 1] - (nfe[j - 1][k - 1] - 1)*p[j - 1] + jitt[j - 1], IloInfinity));
				model.add(IloRange(env, -IloInfinity, sjk[j - 1][k - 1] - alpha[j - 1] - (nfe[j - 1][k - 1] - 1)*p[j - 1] - jitt[j - 1], 0));
				// upper bound
				model.add(IloRange(env, br[j - 1], fe[k - 1] - sjk[j - 1][k - 1], p[j - 1] + jitt[j - 1]));
				
			}
		}


		for (unsigned int k = 1; k <= K; k++) {
			for (unsigned int s = 1; s <= hps.size(); s++) {
				IloExpr S(env); S = sjk[s - 1][k - 1];
				for (unsigned int j = 1; j <= hps.size(); j++) {
					if (s == j) // if s==j then nsjk = nfjk - 1
						model.add(IloRange(env, 0, nsjk[s - 1][j - 1][k - 1] - (nfe[j - 1][k - 1] - 1), 0));
					else {
						//// a coarse upper bound
						model.add(IloRange(env, -IloInfinity, nsjk[s - 1][j - 1][k - 1] - nfe[j - 1][k - 1], 0));

						//// nsjk 
						model.add(IloRange(env, 0, nsjk[s - 1][j - 1][k - 1] * p[j - 1] - (S - alpha[j - 1] - jitt[j - 1]), IloInfinity));
						model.add(IloRange(env, -IloInfinity, nsjk[s - 1][j - 1][k - 1] * p[j - 1] + sigma1 - (S - alpha[j - 1] + jitt[j - 1]) - p[j - 1], 0));
					}
				}
			}
		}

		//// nsb_s_j_k[]

		// Find upperbound
		double ubs = 0;
		for (int i = 0; i < hps.size(); i++) {
			if (p.at(i) + jitt.at(i) > ubs)
				ubs = p.at(i) + jitt.at(i);
		}
		if (pi + ji > ubs)
			ubs = pi + ji;
		/*vector<double> sum_pj = p;
		vector<double> appo_j = jitt;
		appo_j.push_back(ji);
		p.push_back(pi);
		std::transform(sum_pj.begin(), sum_pj.end(), appo_j.begin(), appo_j.end(), std::plus<double>());
		std::sort(sum_pj.begin(), sum_pj.end(),
		[](const double &a, const double &b) { return a > b; });
		double ubs = sum_pj.at(0);*/

		vector< std::vector<std::vector<IloNumVarArray>> > nsb;
		for (unsigned int s = 1; s <= hps.size(); s++) {
			vector < std::vector<IloNumVarArray>> v_;
			for (unsigned int j = 1; j <= hps.size(); j++) {
				vector<IloNumVarArray> v;
				for (unsigned int k = 1; k <= K; k++) {
					IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
					for (unsigned int i = 1; i <= ceil(ubs / p[j - 1]); i++) {
						string name = "nsb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
						vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
					}
					v.push_back(vv);
				}
				v_.push_back(v);
			}
			nsb.push_back(v_);
		}

		////// nfe_j_k = Delta_s + n_s_j_k 
		for (unsigned int s = 1; s <= hps.size(); s++) {
			for (unsigned int j = 1; j <= hps.size(); j++) {
				for (unsigned int k = 1; k <= K; k++) {
					IloExpr Delta(env);
					for (unsigned int i = 1; i <= ceil(ubs / p[j - 1]); i++)
						Delta += nsb[s - 1][j - 1][k - 1][i - 1];
					model.add(IloRange(env, 0, nfe[j - 1][k - 1] - Delta - nsjk[s - 1][j - 1][k - 1], 0));
				}
			}
		}
		////// precedence constraint on ""nsb"s 
		for (unsigned int s = 1; s <= hps.size(); s++) {
			for (unsigned int j = 1; j <= hps.size(); j++) {
				for (unsigned int k = 1; k <= K; k++) {
					for (unsigned int i = 1; i <= ceil(ubs / p[j - 1]); i++) {
						if (i > 1)
							model.add(IloRange(env, 0, nsb[s - 1][j - 1][k - 1][i - 1] - nsb[s - 1][j - 1][k - 1][i - 2], 1));
					}
				}
			}
		}

		//****************************************************//

		//// (5.1.3) nfeb_j_k[] in the interval [a_k-L_k, f_ek]
		vector< std::vector<IloNumVarArray> > nfeb;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			vector<IloNumVarArray> v;
			for (unsigned int k = 1; k <= K; k++) {
				IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
				for (unsigned int i = 1; i <= ceil((Ri + pi + jitt[j - 1]) / p[j - 1]); i++) {
					string name = "nfeb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
					vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
				}
				v.push_back(vv);
			}
			nfeb.push_back(v);
		}
		////// nL_j_k + Delta_e = nfe_j_k 
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K; k++) {
				IloExpr Delta_e(env);
				for (unsigned int i = 1; i <= ceil((Ri + pi + jitt[j - 1]) / p[j - 1]); i++)
					Delta_e += nfeb[j - 1][k - 1][i - 1];
				model.add(IloRange(env, 0, nL[j - 1][k - 1] + Delta_e - nfe[j - 1][k - 1], 0));
			}
		}
		////// precedent constraint on "nfeb"s 
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K; k++) {
				for (unsigned int i = 1; i <= ceil((Ri + pi + jitt[j - 1]) / p[j - 1]); i++) {
					if (i > 1)
						model.add(IloRange(env, 0, nfeb[j - 1][k - 1][i - 1] - nfeb[j - 1][k - 1][i - 2], 1));
				}
			}
		}

		//// nLeb_j_k[] in the interval [f_ek, a_{k+1} - L{k+1}]
		vector< std::vector<IloNumVarArray> > nLeb;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			vector<IloNumVarArray> v;
			for (unsigned int k = 1; k <= K + 1; k++) {
				IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
				for (unsigned int i = 1; i <= ceil((pi - di + ji) / p[j - 1]); i++) {
					string name = "nLeb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
					vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
				}
				v.push_back(vv);
			}
			nLeb.push_back(v);
		}

		////// nfe[k-2]+Delta_p == nL[k-1] 
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 2; k <= K + 1; k++) {
				IloExpr Delta_p(env);
				for (unsigned int i = 1; i <= ceil((pi - di + ji) / p[j - 1]); i++)
					Delta_p += nLeb[j - 1][k - 1][i - 1];
				model.add(IloRange(env, 0, nfe[j - 1][k - 2] + Delta_p - nL[j - 1][k - 1], 0));
			}
		}

		////// nLeb[j][k][i] >= nLeb[j][k][i-1] 
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int k = 1; k <= K + 1; k++) {
				for (unsigned int i = 1; i <= ceil((pi - di + ji) / p[j - 1]); i++) {
					if (i > 1)
						model.add(IloRange(env, 0, nLeb[j - 1][k - 1][i - 1] - nLeb[j - 1][k - 1][i - 2], 1));
				}
			}
		}

		//// nsb_j_k[]
		vector< std::vector<IloNumVarArray> > nsbjk;
		double usb = Ri + pi - bri;
		for (unsigned int j = 1; j <= hps.size(); j++) {
			vector<IloNumVarArray> v;
			for (unsigned int k = 1; k <= K; k++) {
				IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
				for (unsigned int i = 1; i <= ceil(usb / p[j - 1]); i++) {
					string name = "nsb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
					vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
				}
				v.push_back(vv);
			}
			nsbjk.push_back(v);
		}

		////// the cumulative sum of boolean indications
		for (unsigned int j = 1; j <= hps.size(); j++) {
			IloExpr totN(env);
			IloExpr totNe(env);
			for (unsigned int k = 1; k <= K; k++) {

				// enforce f_k
				if (k > 1) {
					for (unsigned int i = 1; i <= ceil((pi) / p[j - 1]); i++)
						totN += nLb[j - 1][k - 1][i - 1];
				}

				model.add(IloRange(env, 0, totN - nL[j - 1][k - 1], 0));

				for (unsigned int i = 1; i <= ceil(ubf / p[j - 1]); i++)
					totN += nfb[j - 1][k - 1][i - 1];

				model.add(IloRange(env, 0, totN - nf[j - 1][k - 1], 0));

				// enforce f_ek
				if (k > 1) {
					for (unsigned int i = 1; i <= ceil((pi - di + ji) / p[j - 1]); i++)
						totNe += nLeb[j - 1][k - 1][i - 1];
				}

				for (unsigned int i = 1; i <= ceil((Ri + pi + jitt[j - 1]) / p[j - 1]); i++)
					totNe += nfeb[j - 1][k - 1][i - 1];
				model.add(IloRange(env, 0, totNe - nfe[j - 1][k - 1], 0));

			}

			// last activation f_k
			for (unsigned int i = 1; i <= ceil((pi) / p[j - 1]); i++)
				totN += nLb[j - 1][K][i - 1];

			model.add(IloRange(env, 0, totN - nL[j - 1][K], 0));

			// last activation f_ek
			for (unsigned int i = 1; i <= ceil((pi - di + ji) / p[j - 1]); i++)
				totN += nLeb[j - 1][K][i - 1];

			model.add(IloRange(env, 0, totNe - nL[j - 1][K], 0));
		}
		////// nLb_j_1_p = 0 [Ex]
		for (unsigned int j = 1; j <= hps.size(); j++) {
			for (unsigned int i = 1; i <= ceil((pi) / p[j - 1]); i++) {
				model.add(IloRange(env, 0, nLb[j - 1][0][i - 1], 0));
			}
		}

		/****** The end of declaration of VARIABLES **********************/




		/********* CONSTRAINTS on schedulability analysis **************/


		// (4.2.1) idle time inside a job window
		for (unsigned int k = 2; k <= K; k++) {
			IloExpr work(env);
			for (unsigned int j = 1; j <= hps.size(); j++) {
				work += (nL[j - 1][k - 1] - nf[j - 1][k - 2]) * c[j - 1];
			}
			model.add(IloRange(env, -IloInfinity, work + idle[k - 2] - (a[k - 1] - L[k - 1] - f[k - 2]), 0));
			model.add(IloRange(env, 0, work + idle[k - 2] - (a[k - 1] - L[k - 1] - f[k - 2]), IloInfinity));
		}


		// (5.1.4) "ft" by accumulating previous workload and idles
		for (unsigned int k = 1; k <= K; k++) {
			IloExpr totI(env);
			for (unsigned int j = 1; j <= hps.size(); j++) {
				totI += nf[j - 1][k - 1] * c[j - 1];
			}

			for (unsigned int h = 1; h <= k; h++)
				totI += delta[h - 1];

			IloExpr totIdle(env);
			if (k > 1)
				for (unsigned int i = 1; i < k; i++) totIdle += idle[i - 1];

			model.add(IloRange(env, 0, totI + totIdle - f[k - 1], 0));
		}

		// (job-kill): "fte" by accummulating previous workload and idles
		for (unsigned int k = 1; k <= K; k++) {
			IloExpr totI(env);
			for (unsigned int j = 1; j <= hps.size(); j++) {
				totI += nfe[j - 1][k - 1] * c[j - 1];
			}

			for (unsigned int h = 1; h <= k; h++) {
				totI += delta[h - 1];
			}

			IloExpr totIdle(env);
			if (k > 1)
				for (unsigned int i = 1; i < k; i++) totIdle += idle[i - 1];

			model.add(IloRange(env, 0, totI + totIdle - fe[k - 1] + b[k - 1] * epsilon, 0));
		}

		// (5.1.5) "ak-Lk" by accumulating previous workload and idles

		for (unsigned int k = 2; k <= K + 1; k++) {
			IloExpr totI(env);
			for (unsigned int j = 1; j <= hps.size(); j++) {
				totI += nL[j - 1][k - 1] * c[j - 1];
			}

			for (unsigned int h = 1; h < k; h++)
				totI += delta[h - 1];

			IloExpr refer(env); refer = a[k - 1] - L[k - 1];

			IloExpr totIdle(env);
			if (k > 1)
				for (unsigned int i = 1; i < k; i++) totIdle += idle[i - 1];

			model.add(IloRange(env, 0, totI + totIdle - refer, 0));
		}

		// (4.2.4) Busy period when bb_{k-1} == 0
		for (unsigned int k = 1; k <= K; k++) {
			IloExpr totI(env);
			for (unsigned int j = 1; j <= hps.size(); j++) {
				totI += (nf[j - 1][k - 1] - nL[j - 1][k - 1]) * c[j - 1];
			}
			totI += delta[k - 1];

			IloExpr refer(env); refer = a[k - 1] - L[k - 1];
			model.add(IloRange(env, -IloInfinity, totI - (f[k - 1] - refer), 0));
			model.add(IloRange(env, 0, totI - (f[k - 1] - refer), IloInfinity));
		}

		// (5.1.6) Busy period (job-kill)
		for (unsigned int k = 1; k <= K; k++) {
			IloExpr totI(env);
			for (unsigned int j = 1; j <= hps.size(); j++) {
				totI += (nfe[j - 1][k - 1] - nL[j - 1][k - 1]) * c[j - 1];
			}
			totI += delta[k - 1] + b[k - 1] * epsilon;

			IloExpr refer(env); refer = a[k -1] - L[k - 1];
			model.add(IloRange(env, -IloInfinity, totI - (fe[k - 1] - refer), 0));
			model.add(IloRange(env, 0, totI - (fe[k - 1] - refer), IloInfinity));

			model.add(IloRange(env, 0, totI - di - sigma1 + (1 - b[k - 1])*M, IloInfinity));
		}

		// ***************** REFINED ANALYSIS ****************************** //
		// (5.2.1) check that sjk is not a valid finish time candidate


		for (unsigned int k = 1; k <= K; k++) {
			for (unsigned int s = 1; s <= hps.size(); s++) {

				IloExpr Sjk(env); Sjk += sjk[s - 1][k - 1];
				IloExpr totI(env);
				for (unsigned int j = 1; j <= hps.size(); j++) {
					totI += nsjk[s - 1][j - 1][k - 1] * c[j - 1];
				}

				for (unsigned int h = 1; h <= k; h++) {
					totI += delta[h - 1];
				}

				IloExpr totIdle(env);
				if (k > 1)
					for (unsigned int i = 1; i < k; i++) totIdle += idle[i - 1];

				model.add(IloRange(env, 0, totI + totIdle - Sjk + b[k - 1] * 1.0 - sigma1, IloInfinity));
			}
		}

		/************************** Additional constraints *************************/
		// C1: constraints on the minimum idle time/workload
		IloExpr idles(env);
		for (unsigned int k = 1; k <= K; k++) {
			idles += idle[k - 1];
			model.add(IloRange(env, idle_lb[k - 1], idles, IloInfinity));
		}

		for (unsigned int j = 1; j <= K; j++) {
			// super coarse upper bound
			model.add(IloRange(env, -IloInfinity, idle[j - 1] - pi + delta[j - 1], 0));
		}

		// C2: let's try to refine the beginning of "bp"
		for (unsigned int k = 2; k <= K + 1; k++) {
			for (unsigned int j = 1; j <= hps.size(); j++)
				model.add(IloRange(env, -IloInfinity, alpha[j - 1] + (nL[j - 1][k - 1] - 1)*p[j - 1] + br[j - 1] - (a[k -1] - L[k - 1]), 0));
		}
		// C3: let's try to refine the beginning of "ft"
		for (unsigned int k = 1; k <= K; k++) {
			for (unsigned int j = 1; j <= hps.size(); j++)
				model.add(IloRange(env, -IloInfinity, alpha[j - 1] + (nf[j - 1][k - 1] - 1)*p[j - 1] + br[j - 1] + sigma1 - f[k - 1], 0));
		}
		// C4: let's try to refine the beginning of "fet"
		for (unsigned int k = 1; k <= K; k++) {
			for (unsigned int j = 1; j <= hps.size(); j++)
				model.add(IloRange(env, -IloInfinity, alpha[j - 1] + (nfe[j - 1][k - 1] - 1)*p[j - 1] + br[j - 1] + sigma1 - fe[k - 1], 0));
		}



		/******* Weakly hard real-time schedulability analysis **********/
		// (4.3) total number of deadline misses
		IloExpr nDmiss(env);
		for (unsigned int k = 1; k <= K; k++)
			nDmiss += b[k - 1];

		model.add(IloRange(env, m + 1, nDmiss, K));
		//model.add(IloRange(env, 0, nDmiss, K));
		model.add(IloMaximize(env, nDmiss));
		/***************************************************************/

		//// Extracting the model
		IloCplex cplex(model);

		const int timeLimit = 60 * 60;
		cplex.setParam(IloCplex::TiLim, timeLimit);

		cplex.setOut(env.getNullStream());

		double ss = cplex.getCplexTime();
		//cplex.exportModel("qcpex1.lp");

		// Stop at the first feasibile solution
		cplex.setParam(IloCplex::IntSolLim, 1);

		// Set maximum number of threads 
		cplex.setParam(IloCplex::Threads, 40);


		cplex.extract(model);
		bool fea = cplex.solve();


		////////////// Analyse data

		const IloAlgorithm::Status solStatus = cplex.getStatus();

		cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
		cout << ">>>>>> " << solStatus << " <<<<<" << std::endl;
		cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;


		///////////////////////////////

		double ee = cplex.getCplexTime();
		bool bounded = false;
		
		swa.time = ee - ss;

		if (fea) {
			swa.bounded = true;
			swa.satisfies_mk = false;
			swa.m = cplex.getObjValue();
		}
		else {
			swa.m = cplex.getObjValue();
			if (swa.time < timeLimit) {
				swa.bounded = true;
				swa.satisfies_mk = true;
			}
			else {
				swa.bounded = false;
				swa.satisfies_mk = false;
			}
		}

		env.end();
		return swa;

	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
		return swa;
	}

	env.end();


}

