/* The state class holds parameters that have been set either on the command line or 
 * in the input file.  
 */
#ifndef GLOBALSTATE_H
#define GLOBALSTATE_H

#define NUM_PARAMS 18
#define FILE_PARAMS 52
#define MAX_CRIT_CATEGORIES 11

#define VL_BINS 9
#define EPI_BINS 8
#define FIRST_BINS 7
#define LAST_BINS 6
#define RISE_BINS 1
#define FALL_BINS 1
#define DAYS_BINS 10
#define PEAK_BINS 10

#define NUM_CRITERIA (VL_BINS+EPI_BINS+FIRST_BINS+LAST_BINS+RISE_BINS+FALL_BINS+DAYS_BINS+PEAK_BINS+2)

#define MAX_X 1000
#define MAX_Y 1000
#define MAX_Z 1000

#define MAX_TIME          10
#define MAX_VIRAL_LOAD    10000000000
#define MAX_INF_CELLS     100000
#define MAX_CD8CELLS      1000000

#define MAX_SAMPLES      20000 /* keep 10 days worth of hex measures */

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "plotpoints.h"
#include "hexcell.h"

typedef struct critType {
	    double low;
	    double high;
	    double mean;
	} critType;

class globalState {
    public:
	gsl_rng * ur;

/* set from input file */
	double TSim;
	double T0_Sim;
	double T2_Sim;
	double T3_Sim;
	double tstep;
	double Tmindel2;
	double Io;
	double To;
	double Yo;
	double Zo;
	double beta_init;
	double beta_low;
	double beta_high;
	double log_p_init;
	double log_p_low;
	double log_p_high;
	double vburstrate_init;
	double vburstrate_low;
	double vburstrate_high;
	double theta_init;
	double theta_low;
	double theta_high;
	double r_init;
	double r_low;
	double r_high;
	double delta_init;
	double delta_low;
	double delta_high;
	double c_init;
	double c_low;
	double c_high;
	double inf_init;
	double inf_low;
	double inf_high;
	double rinf_init;
	double rinf_low;
	double rinf_high;
	double betae_init;
	double betae_low;
	double betae_high;
	double rho_init;
	double rho_low;
	double rho_high;
	double eclipse_init;
	double eclipse_low;
	double eclipse_high;
	double db;
	double N0;
	double an;
	double fpos;
	double swabInterval;
	double statInterval;
	double SnapshotInterval;
	double beta_mean;
	double beta_std;
	double betae_mean;
	double betae_std;
	double vburstrate_mean;
	double vburstrate_std;
	double c_mean;
	double c_std;
	double delta_mean;
	double delta_std;
	double rho_mean;
	double rho_std;
	double theta_mean;
	double theta_std;
	double inf_mean;
	double inf_std;
	double rinf_mean;
	double rinf_std;
	double r_mean;
	double r_std;
	double eclipse_mean;
	double eclipse_std;
	double log_p_mean;
	double log_p_std;

	int PDF_on;
	int gamma_hrs;
	int absorb_hrs;

	/* transmission model variables */
	double log_betaun_init;
	double log_betaun_low;
	double log_betaun_high;

	/* model 6 & 7 vars */
	double gamma_init;
	double gamma_low;
	double gamma_high;
	double gamma_mean;
	double gamma_std;
	double bolus;
	double absorb_init;
	double absorb_low;
	double absorb_high;
	double absorb_mean;
	double absorb_std;
	double Cmax_init;
	double Cmax_low;
	double Cmax_high;
	double Cmax_mean;
	double Cmax_std;
	double IC50_init;
	double IC50_low;
	double IC50_high;
	double IC50_mean;
	double IC50_std;
	double m_init;
	double m_low;
	double m_high;
	double m_mean;
	double m_std;

	double Cmax_0;

	/* model 2 & 3 vars */
	double Iacto;
	double Ilato;
	double alatent;

	/* model 8 vars */
	double kD;
	int cT;
	int nT;

	int N_runs;
	int Total_epis;
	int Tdelay_on;
	int T0_files;
	int use_rho;
	int Tmin;
	int del1;
	int Pulse_neuron;
	int Pulse_regions;
	int Cluster_pulses;
	int Transmission_on;
	int maxCoitalActs;
	int yy;
	int Total_doses;
	int infThreshold;

	int Param_mask;
	int Stop_walk;
	int Bvstop_walk;
	int Max_steps;
	int Printmax;
	int Episode_limit;
	int Sig_test;
	int AutoSnapshot;

	double Tolerance;
	double Size_limit;
	double Input_refresh;

	char *inp_file;
	char *alt_inp_file;

	int Threading;
	int Fit_model;
	int Rand_start;
	int Calc_T0;

	int Model;
	int Model_0;
	int Model_2;
	int Model_3;

	int Regions;
	int Crit_mask;
	int Match_strategy;
	int Search_order; /* 0 - low to high, 1 - start w/ middle then low then high, 2 - high to low */

/* set from criteria file */
	critType crit[NUM_CRITERIA];
	double critWeight[MAX_CRIT_CATEGORIES];

/* used for graphing */
   	int Max_x;
   	int Max_y;
   	int Max_z;

   	int plotVi;
   	int plotVe;
   	int plotInf;
   	int plotCd8;
   	int plotACV;

   	int plotStyle1;
   	int plotStyle2;

   	int plotLogs;
   	int plotColor;
   	int plotRegions;
   	int writeOn;
   	int scrollAxes;

   	int plotOpt1;
   	int plotOpt2;
   	int plotOpt3;
   	int plotOpt4;
   	int plotOpt5;
   	int plotOpt6;
   	int plotOpt7;

 	long int max_vl;
 	long int max_inf;
	long int max_cd8s;
	double max_ACV;

	double time;
	double max_time;
	double plot_bias;
	double plot_span;
	double sampling;
	double refresh;
	double hex_time_bias;
	double crit_start;

	int stopFlag;
	int pauseFlag;
	int Verbose;
	int CritOn;

	FILE *dataF1;
	FILE *dataF2;
	FILE *dataF3;
	FILE *dataF4;
	FILE *dataF5;
	FILE *dataF6;
	FILE *dataF7;
	FILE *dataF8;
	FILE *dataF9;
	FILE *dataF10;
	FILE *dataF11;
	FILE *dataF12;
	FILE *dataF13;

	long time_st;

	plotPoints *points;

 	hexcell **cells;

	int sample_index;
	unsigned long int *vet[MAX_HEXCELLS];
	unsigned long int *vit[MAX_HEXCELLS];
	unsigned long int *inf[MAX_HEXCELLS];
	unsigned long int *cd8[MAX_HEXCELLS];

	double *repro[MAX_HEXCELLS];
	double *diam[MAX_HEXCELLS];
	int *color[MAX_HEXCELLS];

	globalState(void);
	~globalState(void);
};
#endif
