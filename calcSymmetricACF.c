/*	calcSymmetricACF.c
Resvised version of CorrFunc.c for use with Matlab, adapted from:
Correction of the afterpulsing effect in fluorescence correlation spectroscopy using time symmetry analysis
Kunihiko Ishii 1,2 and Tahei Tahara 1,2
1 Molecular Spectroscopy Laboratory, RIKEN, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
2 RIKEN Center for Advanced Photonics, RIKEN, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
2015 Optical Society of America
14 Dec 2015 | Vol. 23, No. 25 | DOI:10.1364/OE.23.032387 | OPTICS EXPRESS 32387
*/


#include <string.h> /* needed for memcpy() */
#include <stdlib.h>
#include "mex.h"


/* mexFunction is the gateway routine for the MEX-file. */
void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

	//  parameters in order found in prhs[]
	unsigned long* taup;			// delay time; here integer number of laser pulses
	unsigned long* dtaup;			// window size; difference between cx[n] and cx[n+1]
	unsigned short* micro;			// microtime
	unsigned long long* macro;		// macrotime
	unsigned short* tcspc_chnsp;	// number of TCSPC channels

	// copy mxArray pointers here for use below - only need to change here
	const mxArray* mxtaup;
	mxtaup = prhs[0];
	const mxArray* mxdtaup;
	mxdtaup = prhs[1];
	const mxArray* mxmicro;
	mxmicro = prhs[2];
	const mxArray* mxmacro;
	mxmacro = prhs[3];
	const mxArray* mxtcspc_chns;
	mxtcspc_chns = prhs[4];

	//  check correct number of dimension of input arrays
	mwSize ndims_micro = mxGetNumberOfDimensions(mxmicro);
	mwSize ndims_macro = mxGetNumberOfDimensions(mxmacro);

	if (ndims_micro != 2 || ndims_macro != 2 ) { 
		mexErrMsgIdAndTxt("calcSymmetricACF:mtime_dims", "Microtime and Macrotime arrays should be Matlab dimensions N x 1, with N photon events.");
	}

	//  check micro and macro have the same length
	const mwSize* microdims;
	const mwSize* macrodims;
	microdims = mxGetDimensions(mxmicro);
	macrodims = mxGetDimensions(mxmacro);
	
	if ( microdims[1] != 1 || macrodims[1] != 1 ) { 
		mexErrMsgIdAndTxt("calcSymmetricACF:mtime_dims", "Microtime and Macrotime arrays should be Matlab dimensions N x 1, with N photon events.");
	}
	if ( microdims[0] != macrodims[0] ) { 
		mexErrMsgIdAndTxt("calcSymmetricACF:mtime_dims", "Microtime and Macrotime arrays should be the same length with N photon events.");
	}
	unsigned long nevents = (unsigned long)microdims[0];

	//		- think about input parameter for tcspc channels, i.e. decay length parameter

	//  TODO - update below prhs indicies follow above changes

	taup =	(unsigned long*)mxGetData(mxtaup);  //  get pointer to tau
	unsigned long tau = *taup;  //  assign tau variable
	dtaup = (unsigned long*)mxGetData(mxdtaup);  //  get pointer to dtau
	unsigned long dtau = *dtaup;  //  assign dtau variable
	tcspc_chnsp = (unsigned short*)mxGetData(mxtcspc_chns);  //  get pointer to tcspc_chns

	micro =	 (unsigned short*)mxGetData(mxmicro);
	macro = (unsigned long long*)mxGetData(mxmacro);

	unsigned long* mp, * mpr, * mp_offset;
	unsigned long i, p, q;
	unsigned short idelay, jdelay;
	unsigned long long t1, t2, t3, t4;

	mwSize m, n;
	m = (mwSize)*tcspc_chnsp;
	n = 1;
	/* create output arrays mp and mp reverse */
	plhs[0] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL); // unsigned long array
	mp = mxGetPr(plhs[0]);
	plhs[1] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL); // unsigned long array
	mpr = mxGetPr(plhs[1]);

	//  Forward
	p = 0;
	q = 0;
	for (i = 0; i < nevents; i++) {
		t1 = macro[i];
		t3 = t1 + tau;
		t4 = t3 + dtau;
		idelay = micro[i];
		mp_offset = mp + idelay;
		while ((p < nevents) && (macro[p] < t3))
			p++;
		while ((q < nevents) && (macro[q] < t4))
			q++;
		*mp_offset += q - p;  // number of events between macrotime Tq and Tp - see equation (3)
	}

	//  Backward
	p = nevents - 1;
	q = nevents - 1;
	for (i = nevents-1; i > 0; i--) {
		t1 = macro[i];
		t3 = t1 - tau;
		t4 = t3 - dtau;
		idelay = micro[i];
		mp_offset = mpr + idelay;
		while ((p > 0) && (macro[p] > t3))
			p--;
		while ((q > 0) && (macro[q] > t4))
			q--;
		*mp_offset += p - q;  // number of events between macrotime Tq and Tp - see equation (3)
	}
	

	return;
}