#ifndef _FTPIC_COMMON_2D_H
#define _FTPIC_COMMON_2D_H

#include <qdsp.h>

extern double DT;
extern double TMAX;

// system size
extern const double XMAX;
extern const double YMAX;

// grid size and spacing
extern int NGRIDX;
extern int NGRIDY;
extern double DX, DY; // spacing

// particle number and properties
extern int PART_NUM;
extern const double PART_MASS;
extern const double PART_CHARGE;
extern const double EPS_0;

extern const double BEAM_SPEED;
extern const double V_TH;

extern double OMEGA_P;

extern const int MODELOG_MAX;
extern FILE *modeLog;
extern int printTime;


int commonInit(int argc, char **argv,
               double **x, double **y, double **vx, double **vy, int **color,
               QDSPplot **xyPlot);

void init2Stream(double *x, double *v, double *vx, double *vy, int *color);

void initLandau(double *x, double *y, double *vx, double *vy, int *color);

void initWave(double *x, double *v, int *color);

#endif
