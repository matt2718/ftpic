#ifndef _QDSP_H
#define _QDSP_H

/* Dummy QDSP header, used for compiling on headless systems without GLFW */

typedef int QDSPplot;

static QDSPplot *qdspInit(const char *title) { return NULL; }
static void qdspDelete(QDSPplot *plot) { return; }
static int qdspUpdate(QDSPplot *plot, double *x, double *y, int *color, int numPoints) { return 0; }
static int qdspUpdateIfReady(QDSPplot *plot, double *x, double *y, int *color, int numPoints) { return 0; }
static int qdspUpdateWait(QDSPplot *plot, double *x, double *y, int *color, int numPoints) { return 0; }
static void qdspRedraw(QDSPplot *plot) { return; }
static void qdspSetFramerate(QDSPplot *plot, double framerate) { return; }
static void qdspSetBounds(QDSPplot *plot, double xMin, double xMax, double yMin, double yMax) { return; }
static void qdspSetConnected(QDSPplot *plot, int connected) { return; }
static void qdspSetPointColor(QDSPplot *plot, int rgb) { return; }
static void qdspSetBGColor(QDSPplot *plot, int rgb) { return; }
static void qdspSetGridX(QDSPplot *plot, double point, double interval, int rgb) { return; }
static void qdspSetGridY(QDSPplot *plot, double point, double interval, int rgb) { return; }

#endif
