#ifndef _mslrdpso_h
#define _mslrdpso_h

#include "gs.h"
#include "ls.h" 
#include "structs.h"

extern FILE *logFile;

// define class for MSLRDPSO algorithmic options
struct MSRD_Options {
	double alpha_start;	// alpha in RDPSO at the beginning of run, generally decreasing linearly
	double alpha_end;	// alpha in RDPSO at the end of run, see mslrdpso.cc
	double beta_start;	// beta in RDPSO at the beginning of run
	double beta_end;	// beta in RDPSO at the end of run
	int FE_interval;      // number of evaluations between two feature exchange operations
	Boole msrd_regenerate_at_limit; // the way to handle particles out of bounds
  // default values for MSRD options :
  public:
    MSRD_Options () :
        alpha_start(0.9),
        alpha_end(0.0),
	beta_start(1.45),
	beta_end(1.00),
        FE_interval(0),
        msrd_regenerate_at_limit(true)
	{ }
	};

class MultiSwarmRDPSO : public Global_Search 
{
	private:
		Population *_Pi;	// best solution for each particle
		Population *_tPi;	// temporary personal best solution used in feature exchange method
		Individual	*_Pg;	// current best solution
		int best; // index in Pi of current global best solution
		unsigned int pop_size;	// population size
		int size;	// number of dimensions (7*nlig + num_torsion)
		float **v;	// velocity
		float *vmax;	//max velocity 
		float *vmin;	// min velocity
		double *xmax;	// max coord bound
		double *xmin;	// min coord bound
		MSRD_Options msrd_options;
	    
		unsigned int generations;
		int outputEveryNgens;
        Output_pop_stats output_pop_stats;
	
		Local_Search *LocalSearchMethod;
	
	public:	
		~MultiSwarmRDPSO();
		MultiSwarmRDPSO(
			float *init_vmax, 
			float *init_vmin, 
			double *init_xmax, 
			double *init_xmin, 
			MSRD_Options init_msrd_options, 
			Local_Search *init_LS, 
			unsigned int init_max_evals, 
			unsigned int init_max_generations, 
			Output_pop_stats output_pop_stats); 			
		
		Individual& getBest();	

      void initialize(const unsigned int, const unsigned int, const int, FILE *);
      unsigned int num_generations(void) const;
		
		// The following part are derived virtual functions
        char *shortname(void);
        char *longname(void);
		void reset(void);
        void reset(const Output_pop_stats&);
        int terminate(void);
        int search(Population &, Eval *, int outlev, FILE * logFile) { return 0; } 
	int localsearch(Population &, Local_Search *, Eval *evaluate, int outlev, FILE * logFile) { return 0; }
        int hysearch(Population &, Local_Search *, int thNum, Eval *, int outlev, FILE * logFile); //MSLDOCK search function
};

inline char * MultiSwarmRDPSO::shortname(void)
{
        return "MSLRDPSO";
}

inline char * MultiSwarmRDPSO::longname(void)
{
        return "Multi-Swarm Lamarckian Random Drift Partile Swarm Optimization";
}

inline MultiSwarmRDPSO::~MultiSwarmRDPSO()
{
	if(_Pi)
		delete _Pi;
	if(_Pg)
		delete _Pg;
	if(v) {
		for(unsigned int i=0;i < pop_size; i++)	
			delete [] v[i];
		delete [] v;
	}
}

inline MultiSwarmRDPSO::MultiSwarmRDPSO(
			float *init_vmax, 
			float *init_vmin, 
			double *init_xmax, 
			double *init_xmin, 
			MSRD_Options init_msrd_options, 
			Local_Search *init_LS, 
			const unsigned int init_max_evals, 
			const unsigned int init_max_generations, 
			Output_pop_stats init_output_pop_stats) :
    Global_Search(  init_max_evals, init_max_generations)
{
	vmax = init_vmax; 
	vmin = init_vmin;			  
	xmax = init_xmax; 
	xmin = init_xmin;			  
	msrd_options = init_msrd_options;
	LocalSearchMethod = init_LS;
    generations = 0;
	output_pop_stats = init_output_pop_stats;
	_Pi =NULL; 
	_Pg = NULL ;
	 v = NULL; 
}

inline Individual& MultiSwarmRDPSO:: getBest()
{
	return *_Pg;
}
inline void MultiSwarmRDPSO::initialize(const unsigned int init_pop_size, const unsigned int ndims, int outlev, FILE *logFile)
{

    pop_size = init_pop_size;
    size = ndims;
}

inline unsigned int MultiSwarmRDPSO::num_generations(void) const
{
   return(generations);
}

// The following part are derived virtual functions
//int search(Population &, Eval *evaluate, int outlev, FILE * logFile);

inline int MultiSwarmRDPSO::terminate(void)
{
   if (max_generations>(unsigned) 0) {
      return((unsigned)generations>=max_generations); 
   } else {
      return(0);  //  Don't terminate
   }
}

	
inline void MultiSwarmRDPSO::reset(void)
{
	generations = 0;
	if(_Pi)
		delete _Pi;
	if(_Pg)
		delete _Pg;
	if(v) {
		for(unsigned int i=0;i < pop_size; i++)	
			delete [] v[i];
		delete [] v;
	}
	
	_Pi = NULL;
	_Pg = NULL;
	v = NULL;
}
	
inline void MultiSwarmRDPSO::reset(const Output_pop_stats &init_output_pop_stats)
{
    output_pop_stats = init_output_pop_stats; 
	reset();
}
	
#define MSRD_D_MAX     (7+MAX_TORS)   // Max number of dimensions of the search space

#endif
