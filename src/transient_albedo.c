#include "transient_albedo.h"

extern double p0;
extern double prev_error;
extern double prev_time;
extern double interror;

const double kp = 0.004;
const double ki = 0.0;
const double kd = 0.001;

int first = 1;

double max(const double a, const double b){ return (a > b) ? a : b; }
double min(const double a, const double b){ return (a < b) ? a : b; }

double transient_albedo(double time, double albedo_coeff, double pboundary)
{
  double alb;
  double dt;
  double derror_dt;
  double error;

  alb = albedo_coeff;

  if (first)
  {
    p0 = pboundary;
    interror = 0.0;
    error = 0.0;
    first = 0;
  }
  else
  {
    error = p0 - pboundary;
    dt = time - prev_time;
    derror_dt = (error - prev_error) / dt;
    interror += error * dt;
    alb += kp * error + kd * derror_dt + ki * interror;
    /* clamp */
    alb = min(max(alb, -1.0), 1.0);
  }

  prev_time = time;
  prev_error = error;

  return alb;
}
