#include "pso.h"
#include "ranlib.h"
#include "support.h" // for evalmode and rep_constants
#include <math.h>
#include <omp.h>

//extern Eval evaluate;
 // MPique - the following static vars should be put in the class...
//static Boole init_links;
//static float prevBestE;
//static	float prevE[PSO_S_MAX];       // E array for particles in previous 
//static	float curE[PSO_S_MAX];        // E array for particles in current step
//static int thNum;
static int cdim;
static int cps;
static int exdim;
static int excps;
static int exf;
static int pexf;
static int movef[PSO_D_MAX];
static int dmat[PSO_D_MAX];
static int cs[PSO_D_MAX];
static int ce[PSO_D_MAX];
static int ds[PSO_D_MAX];
static int de[PSO_D_MAX];

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


/***********************************************************************
 * Particle Swarm Optimization (PSO) with time varying inertia weight w.
 * Global PSO and local S&W search.
 * 
 * Based on PSO in SODOCK: J. Comput Chem, 28, 2, 612-623, 2007.
 * 
 * Huameng Li, 07/2008
 * 
 *
 * Key modifications spring-summer 2011 by R Huey & M Pique at TRSI:
 *   No longer does local search here, the caller will do that
 *     after each generation using our localsearch method, below.
 * Caution: do not reorder the Pop array outside of these functions. MP/RH
 * Packaged algorithm options into a structure to simplify signatures.
 *
 * Caution: not yet thread-safe (local Pi array) and does not delete
 *  the "new"-created arrays, which are members of the ParticleSwarmGS class, see pso.h - M Pique 2014
 ***********************************************************************/
int ParticleSwarmGS::hysearch(Population &Pop, Local_Search *local_method, int thNum, Eval *evaluate, int outlev, FILE * logFile)
{
	//int j;
	//double curVal;
	//double newVal;
	//double Pop_best_value;
	float g_ratio, e_ratio, ratio;
	//int links[PSO_S_MAX][PSO_K_MAX]; // neighbors of i (the [i][0] is always self)
	int cbest[thNum];

	// copy some options for convenience
	float pso_w = pso_options.pso_w;	   // inertia weight
	float w_start = pso_options.pso_w_start;	// pso_w at beginning of run
	float w_end = pso_options.pso_w_end;	// pso_w at conclusion of run
	float c1 = pso_options.c1;
	float c2 = pso_options.c2;
	//int pso_K = pso_options.pso_K;      // number of neighbor particles

	
	// on first call per run, allocate Pi array, allocate and initialize velocity vectors, prevE, curE
	// print pso options if outlev>0
	if(_Pi == NULL) {
#define bool(i)  ((i)?"true":"false")
		
		pr(logFile, "PSO max_generations = %d\n", max_generations);
		pr(logFile, "PSO max_evaluations = %d\n", max_evals);
		pr(logFile, "PSO size = %d\n", size);

		if(outlev>LOGBASIC){
		pr(logFile, "PSO w_start = %.2f\n", pso_options.pso_w_start);
		pr(logFile, "PSO w_end = %.2f\n", pso_options.pso_w_end);
		pr(logFile, "PSO c1 = %.2f\n", pso_options.c1);
		pr(logFile, "PSO c2 = %.2f\n", pso_options.c2);
		pr(logFile, "PSO K = %d\n", pso_options.pso_K);
		pr(logFile, "PSO neighbors_dynamic = %s\n", bool(pso_options.pso_neighbors_dynamic));
		pr(logFile, "PSO neighbors_symmetric = %s\n", bool(pso_options.pso_neighbors_symmetric));
		pr(logFile, "PSO random_by_dimension = %s\n", bool(pso_options.pso_random_by_dimension));
		pr(logFile, "PSO adaptive_velocity = %s\n", bool(pso_options.pso_adaptive_velocity));
		pr(logFile, "PSO regenerate_at_limit = %s\n", bool(pso_options.pso_regenerate_at_limit));
		pr(logFile, "PSO stage2constriction = %s\n", bool(pso_options.pso_stage2constriction));
		pr(logFile, "PSO interpolate_as_scalars = %s\n", bool(pso_options.pso_interpolate_as_scalars));
		}

		fprintf(logFile, "PSO Allocate initial velocity of particles...\n");
		_Pi = new Population(Pop); // copy Pop
		
		pop_size = Pop.num_individuals();
                pr(logFile, "PSO pop_size = %d\n", pop_size); // MP debug

		/*best = 0;
		for(unsigned int i = 1; i < pop_size; i++)
                        if( (*_Pi)[i].value(Normal_Eval) < (*_Pi)[best].value(Normal_Eval) )
                                best = i;
                _Pg = new Individual( (*_Pi)[best] );*/

		// allocate velocity 
                v = new float * [pop_size];
                for(unsigned int i=0; i < pop_size; i++) {
                        v[i] = new float [size];
                }

                // initial velocity of particles
                for(unsigned int i = 0; i < pop_size; i++) {
                        for(int j = 0; j < size; j++) {
                                v[i][j] = random_range(vmin[j], vmax[j]);
                        }
                }
                // MP: note that with adaptive velocity, the BIG will prevent the
                //  first velocity update from occurring
                //for(unsigned int i = 0; i < pop_size; i++) prevE[i] = curE[i] = BIG; // initially unfavorable

                // Display field title
                //pr(logFile, "Generation NumEvals     Pg       PopAvg     PiAvg       w      |V| Avg\n");
                //pr(logFile, "---------- -------- ---------- ---------- ---------- -------- --------\n");                                      
                //pr(logFile, "Generation\tNumEvals\tpi Pop_best\tgBest Pg\tImproved\n");       
                pr(logFile, "Generation NumEvals  gBest E   \n");
                pr(logFile, "---------- -------- ---------- \n");
                //init_links = TRUE;
		
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
	
		exf = 30 * pop_size;
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
	float alpha = 0.0 + (0.9 - 0.0) * (1.0 - ratio);
	float beta = 1.0 + (1.45 - 1.0) * (1.0 - ratio);
	int nexf = pexf + exf;

	int num_tevals[thNum];
	Boole lsf = ranf() < 1/float(thNum);
	pr(logFile, "LS Flag = %s, ", bool(lsf));

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

				//apply local search in the best tPi individual
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
						vmax[j] = (xmax[j] - xmin[j])/2;
						vmin[j] = -vmax[j];
						if(v[i][j] > vmax[j]) v[i][j] = vmax[j];
                        			else if(v[i][j] < vmin[j]) v[i][j] = vmin[j];
					} //next dim
				} //end of velocity update

				// update position (translation, rotation, torsions) Pop[i] of each particle
				for(int i = cs[t]; i <= ce[t]; i++) {
					for(int j = 0; j < size; j++) {
						curVal = Pop[i].phenotyp.gread(j).real;

						// wrap torsions (MP not sure if ever necessary...)
						if(is_conformation_index(j)) curVal=WrpModRad(curVal);

						newVal = curVal + v[i][j];
						
						if(newVal < xmin[j] || newVal > xmax[j]) {
                                			newVal = random_range(xmin[j], xmax[j]);
                        			}
						
						/*if(is_conformation_index(j)) newVal=WrpModRad(newVal);
                        			else if( is_translation_index(j) && (newVal < xmin[j] || newVal > xmax[j]) ) {
                                			if(pso_options.pso_regenerate_at_limit) {
                                				// regenerate at random xyz - do not 
                                				//  change orientation or torsions for now
                                				for(int d=0;d<3;d++) {
                                        				Pop[i].phenotyp.write( random_range(xmin[d], xmax[d]), d);
                                        			}
                                				break;  // end of this individual move
                                			}
                        				else {
                             					if(newVal < xmin[j]) newVal=xmin[j];
                             					else if (newVal > xmax[j]) newVal=xmax[j];
                                			}
						}*/

						// update x,y,z, quaternion, torsion of particle i
						Pop[i].phenotyp.write(newVal, j);
					} // next dim
					// must always normalize Quaternion after modifying its components
					Quat q = Pop[i].phenotyp.readQuat();
					Pop[i].phenotyp.writeQuat( normQuat( q ));
					Pop[i].inverse_mapping(); // MP@@ copy phenotype to genotype (may not be necc)

				}       // end update position

//#pragma omp critical
//	{
				int tbest = cs[t];
                                for(int i = cs[t]+1; i <= ce[t]; i++)
                                        if(Pop[i].value(Normal_Eval) < Pop[tbest].value(Normal_Eval))
                                                tbest = i;

                                //apply local search in the best tPi individual
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
	pr(logFile, "Best value in generation %d, eval %d, value = %f\n",generations, evaluate->evals(), Pg.value(Normal_Eval));

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

		//exchange information between multi swarms
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

	// Output PSO statistics
	if(outputEveryNgens > 0 && 
	  (generations % outputEveryNgens == 0||generations==1)) {
		//pr(logFile, "%d %8d %10.2f %10.2f %10.2f %6.2f %8.2f\n", generations, evaluate->evals(), Pg.value(Normal_Eval), Pop_avg, Pi_avg, w, v_avg);
		//ORIG pr(logFile, "%8d %10ld %10.2f  \n", generations, evaluate->evals(), Pg.value(Normal_Eval));		
		fflush(logFile);
	}
	//pr(logFile, "%8d\t%6d\t%6.3f\t%6.3f\t%6.3f\n", generations, evaluate->evals(), Pop_best_value, Pg.value(Normal_Eval), dist);	
	//pr(logFile, "%6d\t%6d\t%6.3f\t%6.3f\n\n", generations, evaluate->evals(), Pop_best_value, Pg.value(Normal_Eval));		
	//fflush(logFile);
	
	return (0);
}
