#ifndef _FTPIC_COMMON_H
#define _FTPIC_COMMON_H

void init2Stream(double *x, double *v, int *color, int PART_NUM,
                 double XMAX, double BEAM_SPEED);

void initDisplace(double *x, double *v, int *color, int PART_NUM,
                  double XMAX, double OMEGA_P);

void initRhoSin(double *x, double *v, int *color, int PART_NUM, double XMAX);

#endif
