/*=================================================================
 * BaseGenomeBreaks.h
 *
 *=================================================================*/
/*
	This File is part of GADA

	GADA v1.0 Genome Alteration Detection Algorithm
    Copyright (C) 2008  Childrens Hospital of Los Angeles

	GADA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    GADA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GADA.  If not, see <http://www.gnu.org/licenses/>.

	Author:
		Roger Pique-Regi    piquereg@usc.edu

*/

#ifndef _BaseGADA_H_
#define _BaseGADA_H_

//#include "matlabdefines.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#ifndef GADABIN
#include <boost/python.hpp>
#endif
using namespace std;


void reconstruct (double *wr,long M,double *aux_vec);
void BubbleSort (long *I,long L);
void doubleBubbleSort (double *D,long *I,long L);
void TrisolveREG(double *t0,double *tu,double *tl,double *coef,double *sol,long sizeh0);
void DiagOfTriXTri(double *ll,double *l0,double *lu,double *rl,double *r0,double *ru,double *d,long N);
void tridiagofinverse(double *t0,double *tl,double *itl,double *it0,double *itu,long N,double *d,double *e);
void ForwardElimination(double *A,long N);
void BackSubstitution(double *A,long N);
void BackwardElimination(double *A,long N);
void TriSolveINV(double *AA,long M, long N, double *x,double *d,double *e);
void ComputeH(double *h0, double *h1, long M);
void ComputeFdualXb(long M, double *b);
// 20080119 REMOVED void ComputeHs(long *s,double *a,long M,long Ms,double *h0,double *h1);
void ComputeHs(long *s,long M,long Ms,double *h0,double *h1);
void TriSymGaxpy(double *t0, double *t1, double *x, long M, double *y);
void ComputeT(double *h0,double *h1,long M,double *alfa,double sigma,double *t0,double *tl,double *tu);
long findminus(double *alpha,long Ms,double maxalpha,long *sel);
long simpletresholding(double *inputvector,long N,double thres,double *disc);
void computesegmentmeans(double *inputvector,long N,double *disc,long numdisc,double *amp);
void reconstructoutput(double *rec,long N,double *disc,long numdisc,double *amp);
long SBL(
    double *y, //I -- 1D array with the input signal
    long *I, //IO -- 1D array with the initial (final) candidate breakpoints
    double *alpha, //I -- 1D array with the initial (final) hyperparameter inv. varainces.
    double *w, //O -- 1D array containing the breakpoint weigths or posterior mean.
    double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha and H
    long M, //Initial size of the array in y
    long *K, //Size of the I alpha w

    //Algorithm parameters:
    double sigma2, //Noise estimated
    double a,      //
    double b,
    double maxalpha,  //Basis reduction parameter
    long    maxit,     //Max number of iterations
    double tol,       //Tolerance for convergence
    long debug       //verbosity... set equal to 1 to see messages  0 to not see them
    );

long BEthresh( //To eliminate...
    double *Scores,
    long Nscores,
    double *wr,
    long *indsel,
    long *pointNumRem,
    double *pointTau
    );

long SBLandBE( //Returns breakpoint list lenght.
    double *tn,
    long M,  //length of the noisy signal tn
    double *sigma2, //If sigma2 < 0, compute sigma2 (Input/Output)
    double a,      // SBL parameter
    double T,      // Threshold to prune
    long MinSegLen,	//Minimum length of the segment.
    long **pI,		//Returns breakpoint positions
    double **pw,	//Returns breakpoint weights.
    long debug		//verbosity... set equal to 1 to see messages  0 to not see them
    //long *pK
    );


void Project(
    double *y,
    long M,
    long *I,
    long L,
    double *xI,
    double *wI
    );
void IextToSegLen(
	long *Iext, // Enters Iext
	long *SegLen,		// Outputs SegLen.. can be the same as Iext?
	long K			// Length Iext - 1
	);
void IextWextToSegAmp(
	long *Iext,
	double *Wext,
	double *SegAmp,
	long K
	);
void CompZ(// computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
	double *y,
	double *z,
	long M
	);
void ComputeHsIext(
    //input variables:
    long *Iext,     // Indices of selection,
    long K,     // Length of the indices,
    double *h0, // Returning diagonal of H,
    double *h1  // Returning upper diagonal of H
    );
void ProjectCoeff ( //IextYobs2Wext
    double *y,
    long M,
    long *Iext,
    long K,
    double *Wext
	);
void
CollapseAmpTtest(//Uses a T test to decide which segments collapse to neutral
	double *SegAmp, //Segment Amplitudes (input output)
	const long *SegLen, //Segment Lengths
	long K, //Number of segments.
	double BaseAmp, //Reference amplitude to compare
	double sigma2,  //Reference noise
	double T		//Critical value that decides when to colapse
	);
double // Returns BaseAmp corresponding to the base level.
CompBaseAmpMedianMethod( //Computes the median recontruction level, as baseline level.
	const long *SegLen,    //Lengths corresponding to the amplitudes
	const double *SegAmp, //Amplitudes !!! assumed already ordered...
	long K
	);

void
ClassifySegments(
	double *SegAmp,
	long *SegLen,
	double *SegState,
	long K,
	double BaseAmp,
	double ploidy,
    double sigma2,  //Reference noise
	double T		//Critical value that decides when to colapse
	);

void
ComputeTScores(
	const double *Wext,
	const long *Iext,
	double *Scores,
	long K,
	long start,
	long end
	);

long BEwTscore(
    double *Wext,  //IO Breakpoint weights extended notation...
	long *Iext,     //IO Breakpoint positions in extended notation...
	double *tscore,
    long *pK,       //IO Number breakpoint positions remaining.
    double T      //IP  Threshold to prune
	);

long BEwTandMinLen( //Returns breakpoint list lenght. with T and MinSegLen
    double *Wext,  //IO Breakpoint weights extended notation...
	long *Iext,     //IO Breakpoint positions in extended notation...
    long *pK,       //IO Number breakpoint positions remaining.
	double sigma2, //IP If sigma2
    double T,      //IP  Threshold to prune,  T=T*sqrt(sigma2);
    long MinSegLen  //IP Minimum length of the segment.
);
long
RemoveBreakpoint(
	double *Wext,
	long *Iext,
	long K,
	long jrem
	);


#endif //_BaseGADA_H_
