#include "mslrdpso.h"
#include "ranlib.h"
#include "support.h" // for evalmode and rep_constants
#include <math.h>
#include <omp.h>

/***************************************************************************************************
 * Random Drift Particle Swarm Optimization (RDPSO) with a multi-swarm strategy,
 * and a modified implementation of S&W local search method to adapt to the  multi-swarm structure
 * 
 * Based on the paper by Chao Li, Jun Sun, Xiaojun Wu, and Vasile Palade:
 * MSLDOCK: multi-swarm optimization algorithm for flexible ligand docking 
 *
 * The code is written by
 * Chao Li, 02/2020
 * 
 * This method is thread-safe and thus can be implemented in parallel by using OpenMP
 * Set msldock_run in .dpf file to call this function
 * Caution: this version of MSLDOCK only supports the command-line usage
 * Running in serial: msldock -p file.dpf -l file.dlg
 * Running in parallel: msldock.omp -p file.dpf -l file.dlg 
 * (in parallel mode, number of threads equals to the number of sub-swarms if possible)
 **************************************************************************************************/

static int cdim; //the feature dimension each sub-swarm contains
static int cps; //the particles each sub-swarm contains 
static int exdim; //the remaining feature dimension shared by the first sub-swarms equally
static int excps; //the remaining particles shared by the first sub-swarms equally
static int exf; //number of evaluations between two feature exhange method
static int pexf; //number of evaluations of the pervious feature exhange
static int movef[MSRD_D_MAX]; //identify whether move features or not
static int dmat[MSRD_D_MAX]; //identify each feature belongs to which sub-swarm
static int cs[MSRD_D_MAX]; //start sequence number of particles for each sub-swarm
static int ce[MSRD_D_MAX]; //end sequence number of particles for each sub-swarm
static int ds[MSRD_D_MAX]; //start sequence number of feature dimensions for each sub-swarm
static int de[MSRD_D_MAX]; //end sequence number of feature dimensions for each sub-swarm

inline float Norm(float *x, int n)
{
	float s = 0;
	for(int i=0;i<n;i++)
		s += (x[i]*x[i]);
	return sqrt(s);
}

void initArray(int *array, int len)
{
	for(int i = 0; i < len; i++) 
		array[i] = i;
}

int randin(int low, int high)
{
	return floor(low+(high-low)*ranf());
}


int MultiSwarmRDPSO::hysearch(Population &Pop, Local_Search *local_method, int thNum, Eval *evaluate, int outlev, FILE * logFile)
{
	float g_ratio, e_ratio, ratio;

	// copy some options for convenience
	float alpha_start = msrd_options.alpha_start;	// alpha in RDPSO at the beginning of run, generally decreasing linearly
	float alpha_end = msrd_options.alpha_end;	// alpha in RDPSO at the end of run
	float beta_start = msrd_options.beta_start;	// beta in RDPSO at the beginning of run
	float beta_end = msrd_options.beta_end;		// beta in RDPSO at the end of run
	int exf = msrd_options.FE_interval;      // number of evaluations between two feature exchange operations

	if(exf == 0) {
		exf = 30 * pop_size;
	}
	
	int cbest[thNum]; //record the sequence number of the best Pi for each sub-swarm

	// on first call per run, allocate Pi array
	// print MSLDOCK options if outlev>0
	if(_Pi == NULL) {
#define bool(i)  ((i)?"true":"false")
		
		pr(logFile, "MSLDOCK max_generations = %d\n", max_generations);
		pr(logFile, "MSLDOCK max_evaluations = %d\n", max_evals); //generally use max_evaluations as termination condition
		pr(logFile, "MSLDOCK size = %d\n", size);

		pr(logFile, "MSLDOCK number of sub-swarms = %d\n", thNum);
		pr(logFile, "MSLDOCK FE_interval = %d (default:30*pop_size)\n", exf);

		if(outlev>LOGBASIC){
		pr(logFile, "MSLDOCK alpha_start = %.2f\n", msrd_options.alpha_start);
		pr(logFile, "MSLDOCK alpha_end = %.2f\n", msrd_options.alpha_end);
		pr(logFile, "MSLDOCK beta_start = %.2f\n", msrd_options.beta_start);
		pr(logFile, "MSLDOCK beta_end = %.2f\n", msrd_options.beta_end);
		pr(logFile, "MSLDOCK regenerate_at_limit = %s\n", bool(msrd_options.msrd_regenerate_at_limit));
		}

		_Pi = new Population(Pop); // copy Pop
		
		pop_size = Pop.num_individuals();
                pr(logFile, "MSLDOCK pop_size = %d\n", pop_size); // MP debug

		// allocate velocity 
                v = new float * [pop_size];
                for(unsigned int i=0; i < pop_size; i++) {
                        v[i] = new float [size];
                }

                // initial velocity of particles
                //for(unsigned int i = 0; i < pop_size; i++) {
                //        for(int j = 0; j < size; j++) {
                //                v[i][j] = random_range(vmin[j], vmax[j]);
                //        }
                //}
                //pr(logFile, "Generation NumEvals  gBest E   \n");
                //pr(logFile, "---------- -------- ---------- \n");
		
		//thNum = omp_get_max_threads();
                //if (thNum > size) thNum = size;
		//if (thNum > floor(size/2)) thNum = floor(size/2);
		//if (thNum < 2) thNum = 2;
                cdim = floor(size/thNum);
                exdim = size % thNum;
                cps = floor(pop_size/thNum);
                excps = pop_size % thNum;
		_tPi = new Population(Pop);

		initArray(dmat, size);
		genprm(dmat, size);

#pragma omp parallel
	{
#pragma omp for
		for(int t = 0; t < thNum; t++) {
			if(t < exdim) {
				ds[t] = t * (cdim+1);
				de[t] = (t+1) * (cdim+1) - 1;
			}
			else {
				ds[t] = t * cdim + exdim;
				de[t] = (t+1) * cdim + exdim - 1;
			}
			if(t < excps) {
				cs[t] = t * (cps+1);
				ce[t] = (t+1) * (cps+1) - 1;
			}	
			else {
				cs[t] = t * cps + excps;
				ce[t] = (t+1) * cps + excps - 1;
			}
			cbest[t] = cs[t];
			//find best Pi for each sub-swarm
			for(int i = cs[t]+1; i <= ce[t]; i++)
                        	if( (*_Pi)[i].value(Normal_Eval) < (*_Pi)[cbest[t]].value(Normal_Eval) )
                                	cbest[t] = i;
		}
	}
//end of parallel
	
		int tbest = 0;
		for(int t = 1; t < thNum; t++)
                	if((*_Pi)[cbest[t]].value(Normal_Eval) < (*_Pi)[cbest[tbest]].value(Normal_Eval) )
               			tbest = t;
		best = cbest[tbest];
        	_Pg = new Individual( (*_Pi)[best] );

		pexf = 0;
		for(int t = 0; t < thNum; t++)
			movef[t] = 0;
	}

	Population &Pi = (Population &)(*_Pi);
	Population &tPi = (Population &)(*_tPi);
	Individual &Pg = (Individual &)(*_Pg);

	g_ratio = max_generations>0 ? (float)generations / max_generations : 0;
	e_ratio = max_evals>0 ? (float)evaluate->evals() / max_evals : 0;
	ratio =  max(g_ratio, e_ratio);
	float alpha = alpha_end + (alpha_start - alpha_end) * (1.0 - ratio);
	float beta = beta_end + (beta_start - beta_end) * (1.0 - ratio);
	int nexf = pexf + exf;

	int num_tevals[thNum];
	Boole lsf = ranf() < 1/float(thNum); //in each generation has a probability to apply local search to all sub-swarms

#pragma omp parallel
        {
#pragma omp for
		for(int t = 0; t < thNum; t++) {
			if(movef[t] == 1) {
//#pragma omp critical
//	{
				Eval tmpEval = *evaluate;
				for(int i = cs[t]; i <= ce[t]; i++) {
					tPi[i].phenotyp.pevaluate = &tmpEval;
				}

				int tbest = cs[t];
				for(int i = cs[t]+1; i <= ce[t]; i++)
					if(tPi[i].value(Normal_Eval) < tPi[tbest].value(Normal_Eval))
						tbest = i;

				//apply local search in the best tPi in each sub-swarm
				if(local_method != NULL && lsf) {
					local_method->search(tPi[tbest], &tmpEval, outlev, logFile);
				}
//	}

				for(int i = cs[t]; i <= ce[t]; i++ ) {
                			if( tPi[i].value(Normal_Eval) < Pi[i].value(Normal_Eval) )
                        			Pi[i] = tPi[i];
				}
        			cbest[t] = cs[t];
        			for(int i = cs[t]+1; i <= ce[t]; i++ ) 
                			if( Pi[i].value(Normal_Eval) < Pi[cbest[t]].value(Normal_Eval) ) 
                        			cbest[t] = i;
				movef[t] = 0;

				num_tevals[t] = tmpEval.evals() - evaluate->evals();
			}
			else {
				Eval tmpEval = *evaluate;
				for(int i = cs[t]; i <= ce[t]; i++) {
					Pop[i].phenotyp.pevaluate = &tmpEval;
				}

				// update velocity
				double curVal, newVal;
				double Pavg[size];
				for(int j = 0; j < size; j++) {
					Pavg[j] = 0.0;
					for(int i = cs[t]; i <= ce[t]; i++) {
						Pavg[j] += Pi[i].phenotyp.gread(j).real;
					}
					Pavg[j] /= (ce[t] - cs[t] + 1);
				}
				for(int i = cs[t]; i <= ce[t]; i++) {
					double fi, ri, pp;
					for(int j = 0; j < size; j++) {
						fi = ranf();
						ri = gennor(0,1);
						curVal = Pop[i].phenotyp.gread(j).real;
						pp = fi * Pi[i].phenotyp.gread(j).real
						   + (1-fi) * Pi[cbest[t]].phenotyp.gread(j).real;
						v[i][j] = alpha * ri * fabs(Pavg[j] - curVal) + beta * (pp - curVal);
						//vmax[j] = (xmax[j] - xmin[j])/2;
						//vmin[j] = -vmax[j];
						if(v[i][j] > vmax[j]) v[i][j] = vmax[j];
                        			else if(v[i][j] < vmin[j]) v[i][j] = vmin[j];
					} //next dim
				} //end of velocity update

				// update position (translation, rotation, torsions) Pop[i] of each particle
				for(int i = cs[t]; i <= ce[t]; i++) {
					for(int j = 0; j < size; j++) {
						curVal = Pop[i].phenotyp.gread(j).real;

						// wrap torsions (not sure if ever necessary...)
						if(is_conformation_index(j)) curVal=WrpModRad(curVal);

						newVal = curVal + v[i][j];
						
						/*if(newVal < xmin[j] || newVal > xmax[j]) {
                                			newVal = random_range(xmin[j], xmax[j]);
                        			}*/
						
						if(is_conformation_index(j)) newVal=WrpModRad(newVal);
                        			else if( is_translation_index(j) && (newVal < xmin[j] || newVal > xmax[j]) ) {
                                			if(msrd_options.msrd_regenerate_at_limit) {
                                				//regenerate at random xyz - do not change orientation or torsions
                                				for(int d=0;d<3;d++) {
                                        				Pop[i].phenotyp.write( random_range(xmin[d], xmax[d]), d);
                                        			}
                                				break;  // end of this particle move
                                			}
                        				else {
                             					if(newVal < xmin[j]) newVal=xmin[j];
                             					else if (newVal > xmax[j]) newVal=xmax[j];
                                			}
						}

						// update x,y,z, quaternion, torsion of particle i
						Pop[i].phenotyp.write(newVal, j);
					} // next dim
					// must always normalize Quaternion after modifying its components
					Quat q = Pop[i].phenotyp.readQuat();
					Pop[i].phenotyp.writeQuat( normQuat( q ));
					Pop[i].inverse_mapping(); // copy phenotype to genotype (may not be necc)

				}       // end update position

//#pragma omp critical
//	{
				int tbest = cs[t];
                                for(int i = cs[t]+1; i <= ce[t]; i++)
                                        if(Pop[i].value(Normal_Eval) < Pop[tbest].value(Normal_Eval))
                                                tbest = i;

                                //apply local search in the best Pop in each sub-swarm
                                if(local_method != NULL && lsf) {
                                	local_method->search(Pop[tbest], &tmpEval, outlev, logFile);
                                }
//	}
				//update Pi, cbest
				for(int i = cs[t]; i <= ce[t]; i++ ) {
                                        if( Pop[i].value(Normal_Eval) < Pi[i].value(Normal_Eval) )
                                                Pi[i] = Pop[i];
                                }
                                cbest[t] = cs[t];
                                for(int i = cs[t]+1; i <= ce[t]; i++ )
                                        if( Pi[i].value(Normal_Eval) < Pi[cbest[t]].value(Normal_Eval) )
                                                cbest[t] = i;

				num_tevals[t] = tmpEval.evals() - evaluate->evals();
			}
                }
        }
//end of parallel

	for(int t = 0; t < thNum; t++)
		evaluate->num_evals += num_tevals[t];

	//find best and Pg from Pi[cbest]
	int tbest = 0;
	for(int t = 1; t < thNum; t++)
		if(Pi[cbest[t]].value(Normal_Eval) < Pi[cbest[tbest]].value(Normal_Eval) )
			tbest = t;
        best = cbest[tbest];
	Pg = Pi[best];
	if(outlev>LOGRUNV) {
		pr(logFile, "In generation %d, LS_flag = %s, eval %d, best energy value = %f\n",generations, bool(lsf), evaluate->evals(), Pg.value(Normal_Eval));
	}
	
	//shuffle dmat & exchange informtion
	if(evaluate->evals() >= nexf) {
		pexf = nexf;
		for(int t = 0; t < thNum; t++) 
			movef[t] = 1;

		//random shuffle dmat
		int ndmat[size]; //dmat of next generation, an array
		int prevc[size]; //record the dimension in which multi swarm
		int rindex[size];
		initArray(rindex, size);
		genprm(rindex, size);
		for(int j = 0; j < size; j++) {
			ndmat[j] = dmat[rindex[j]];
			for(int t = 0; t < thNum; t++) {
				if( rindex[j] >= ds[t] && rindex[j] <= de[t])
					prevc[j] = t;
			}
		}

		//feature exchange method applied between sub-swarms
		for(int i = 0; i < pop_size; i++) 
			tPi[i] = Pi[i];

#pragma omp parallel
	{
#pragma omp for
		for(int t = 0; t < thNum; t++){
			for(int j = ds[t]; j <= de[t]; j++) {
				double exVal;
				dmat[j] = ndmat[j];
				int dcp = (prevc[j] < excps) ^ (t < excps);
				if( dcp == 0) {
					int tline;
					if(t < excps) 
						tline = cps + 1;
					else
						tline = cps;
					int rindex1[tline];
					initArray(rindex1, tline);
					genprm(rindex1, tline);
					for(int i = 0; i < tline; i++) {
						exVal = Pi[ cs[prevc[j]] + rindex1[i] ].phenotyp.gread(dmat[j]).real;
						tPi[ cs[t] + i ].phenotyp.write(exVal, dmat[j]);
					}
					initArray(rindex1, tline);
					genprm(rindex1, tline);
					for(int i = 0; i < tline; i++) {
						exVal = Pi[ cs[t] + rindex1[i] ].phenotyp.gread(dmat[j]).real;
						tPi[ cs[prevc[j]] + i ].phenotyp.write(exVal, dmat[j]);
					}
				}
				else {
					if(prevc[j] > t) {
						int cpyline = randin(cs[prevc[j]], ce[prevc[j]]+1);
						int rindex1[cps];
						initArray(rindex1, cps);
						genprm(rindex1, cps);
						for(int i = 0; i < cps; i++) {
							exVal = Pi[ cs[prevc[j]] + rindex1[i] ].phenotyp.gread(dmat[j]).real;
							tPi[ cs[t] + i ].phenotyp.write(exVal, dmat[j]);
						}
						exVal = Pi[cpyline].phenotyp.gread(dmat[j]).real;
						tPi[ce[t]].phenotyp.write(exVal, dmat[j]);
						int rindex2[cps+1];
						initArray(rindex2, cps+1);
						genprm(rindex2, cps+1);
						for(int i = 0; i < cps; i++) {
							exVal = Pi[ cs[t] + rindex2[i] ].phenotyp.gread(dmat[j]).real;
							tPi[ cs[prevc[j]] + i ].phenotyp.write(exVal, dmat[j]);
						}
					}
					else {
						int rindex2[cps+1];
						initArray(rindex2, cps+1);
						genprm(rindex2, cps+1);
						for(int i = 0; i < cps; i++) {
							exVal = Pi[ cs[prevc[j]] + rindex2[i] ].phenotyp.gread(dmat[j]).real;
							tPi[ cs[t] + i ].phenotyp.write(exVal, dmat[j]);
						}
						int cpyline = randin(cs[t], ce[t]+1);
						int rindex1[cps];
						initArray(rindex1, cps);
						genprm(rindex1, cps);
						for(int i = 0; i < cps; i++) {
							exVal = Pi[ cs[t] + rindex1[i] ].phenotyp.gread(dmat[j]).real;
							tPi[ cs[prevc[j]] + i ].phenotyp.write(exVal, dmat[j]);
						}
						exVal = Pi[cpyline].phenotyp.gread(dmat[j]).real;
						tPi[ce[prevc[j]]].phenotyp.write(exVal, dmat[j]);
					}
				}
			}
			for(int i = cs[t]; i <= ce[t]; i++) {
				Quat q = tPi[i].phenotyp.readQuat(); 
				tPi[i].phenotyp.writeQuat( normQuat( q ));
				tPi[i].inverse_mapping();
			}
		}
	}
//end of parallel

	}
	
	if(evaluate->evals() >= max_evals || generations >= max_generations) {
 		for(unsigned int i = 0; i < pop_size; i++) {
 			Pop[i]=Pi[i];
 		}
 	}

	generations++;

	if(outputEveryNgens > 0 && 
	  (generations % outputEveryNgens == 0||generations==1)) {
		fflush(logFile);
	}
	
	return (0);
}
