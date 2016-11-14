#pragma once
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
//#include <tchar.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
//#include "parameter.h"

#define NR_END 1
#define FREE_ARG char*
#define MAXITS 1000
#define EPS 1.0e-5

void nrerror(char error_text[]);

void tridag(double a[], double b[], double c[], double r[], double u[],
	unsigned long n);

void formdx(double *dx);

void formdr(double *dr);

void formABA1(double *A, double *B, double *A1, double dx[],double Dc);

void formCDC1(double *C, double *D, double *C1, double dr[],double Dc);

double FdCa(double Cc, double Ca, double ka, double ka1, double kd);

double FKa(double Ca,double kc,double Kca);

double FKm(double Mt, double kp1, double KMT);
