#ifndef TRANSIENT_ALBEDO_H
#define TRANSIENT_ALBEDO_H

extern double p0;
extern double prev_error;
extern double prev_time;
extern double interror;

extern const double kp;
extern const double ki;
extern const double kd;

extern int first;

double max(const double a, const double b);
double min(const double a, const double b);

double transient_albedo(const double time, const double albedo_coeff, const double pboundary);

#endif
