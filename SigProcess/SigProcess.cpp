// SigProcess.cpp : Defines the entry point for the console application.
// Written by Jaafar Alsalaet
//Last edited on 30-8-2021
//c++ implementation of spectral correlation/coherence
//using ACP, FastACP, FAM, FSC and fast Dirichlet kernel SC
//
//FFTW library is required to compile and run the project. Depending on your buil configuration CPU type, 
//rename either libfftw3-3-x86.dll or libfftw3-3-x64.dll to libfftw3-3.dll and put it in your application .exe folder

#include "stdafx.h"
#include "sigpro.h"
#include <fftw3.h>
#include <windows.h>
#include <sstream>
#define DbgMsg( s )            \
{                             \
std::wostringstream os_;    \
   os_ << s;                   \
   OutputDebugStringW( os_.str().c_str() );  \
}

const double pi = 3.1415926535897932384626433832795;

polar operator +(polar const& c1, polar const& c2) {
	double real, imag;
	real = c1.amp * cos(c1.angle) + c2.amp * cos(c2.angle);
	imag = c1.amp * sin(c1.angle) + c2.amp * sin(c2.angle);	
	return polar(sqrt(imag*imag + real*real), atan2(imag, real));;
}

polar operator -(polar const& c1, polar const& c2) {
	complex c3;
	c3.real = c1.amp * cos(c1.angle) - c2.amp * cos(c2.angle);
	c3.imag = c1.amp * sin(c1.angle) - c2.amp * sin(c2.angle);
	return c3.toPolar();
}

polar operator *(polar const& c1, polar const& c2) {
	polar c3;
	c3.amp = c1.amp * c2.amp;
	c3.angle = c1.angle + c2.angle;
	return c3;
}

polar operator /(polar const& c1, polar const& c2) {
	polar c3;
	c3.amp = c1.amp / c2.amp;
	c3.angle = c1.angle - c2.angle;
	return c3;
}

complex operator +(complex const& c1, complex const& c2) {
	complex c3;
	c3.real = c1.real + c2.real;
	c3.imag = c1.imag + c2.imag;
	return c3;
}

complex operator -(complex const& c1, complex const& c2) {
	complex c3;
	c3.real = c1.real - c2.real;
	c3.imag = c1.imag - c2.imag;
	return c3;
}

complex operator *(complex const& c1, complex const& c2) {
	complex c3;
	c3.real = c1.real*c2.real - c1.imag*c2.imag;
	c3.imag = c1.real*c2.imag + c2.real*c1.imag;
	return c3;
	//return complex(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c2.real*c1.imag);
}
complex operator *(float const& c1, complex const& c2) {
	complex c3;
	c3.real = c1 * c2.real;
	c3.imag = c1 * c2.imag;
	return c3;	
}
complex operator *(double const& c1, complex const& c2) {
	complex c3;
	c3.real = c1 * c2.real;
	c3.imag = c1 * c2.imag;
	return c3;
}
complex operator /(complex const& c1, complex const& c2) {
	complex c3;
	double r2;
	r2 = c2.real*c2.real + c2.imag*c2.imag;
	c3.real = (c1.real*c2.real + c1.imag*c2.imag) / r2;
	c3.imag = (c2.real*c1.imag - c1.real*c2.imag) / r2;
	return c3;
}

void fft(complex *A, int Nb, int stat)
{
	//stat=1 forward transform, -1 backward transform
	// starting index for the input and output data stream is 0
	int lpk, l, M, me1, k, j, nbd2, nbm1, N;
	complex u1, w1;
	complex to1;
	//float uu2 ;

	N = (int)(log2((double)Nb));

	//		BitReversal(A, N); //~same speed as the following code

	nbd2 = Nb / 2;
	nbm1 = Nb - 1;
	j = 0;	
	for (l = 0; l < nbm1; l++)
	{
		if (l < j)
		{
			to1 = A[j];
			A[j] = A[l];
			A[l] = to1;
		}

		k = nbd2;
	label3:
		if (k > j)
		{
			j = j + k;
		}
		else
		{
			j = j - k;
			k = k / 2;
			goto label3;
		}
	}

	for (M = 1; M <= N; M++)
	{
		u1 = complex(1.0, 0.0);
		me1 = 1 << M;
		k = me1 >> 1;
		w1 = complex(cos(pi / k), -1.0 * stat * sin(pi / k));
		for (j = 0; j < k; j++)
		{
			for (l = j; l< Nb; l += me1)
			{
				lpk = l + k;
				to1 = A[lpk] * u1;
				A[lpk] = A[l] - to1;
				A[l] = A[l] + to1;

			}//end of L			
			u1 = u1 * w1;
		} //end of j
	} //end {of m}

}
void ffthelper(complex* A, complex* B, int N)
{//helper function 
	int i;
	for (i = 0;i < N;i++)
		B[i] = A[i];
	fft(B, N, 1);
}
int FSCoh(float* x, float &da, int a1, int a2, int L, int Noverlap, int Nw, int opt, float* SCo)
{


	/*
	calculate spectral correlation/coherance using Fast ACP

	da: delta alpha relative to Fs
	x() input data
	a1: minimum cyclic freq, must be 0, this is ignored input just to keep compatibility with other functions form
	a2: maximum cyclic frequency count such that max. alpha = da * a2
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data: of size Nw/2 * a2
	*/
	int nfft, hnfft, KN, i, k, j, R, Nb, ri, ci, mri, inc;
	bool Cohr;
	double mk, cp1;
	R = Nw - Noverlap;
	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out;
	fftw_plan p;
	complex *sptr;
	//auto a = new complex[Nw];
	//auto b = new complex[Nw];
	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	p = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw;i++)
		windat[i] = 2.0 * (0.5 - 0.5 * cos(2 * i * pi / Nw));

	KN = (int)((L - Noverlap) / (double)(Nw - Noverlap));
	Cohr = false;
	if ((opt == 1) || (opt == 3)) Cohr = true;
	//std::cout << "KN: " << KN << "\n";	
	auto dshift = new complex[Nw];
	complex sfa, sfm, kfk;
	auto CPY2 = new double[hnfft](); //putting () will initialize all elements to zero
	auto SCoCmplx = new complex[hnfft * a2];
	auto SCoX2 = new double[hnfft * a2]();
	double KMU = 0.0;
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = KN * KMU;
	Nb = (int)round((1.0f / da / Nw));
	DbgMsg("FastACP Nb: " << Nb << "\n");
	if (Nb < 1) {
		//("Window length must be smaller than signal length");
		return -1;
	}
	da = 1.0F / (Nw * Nb);
	auto Kft = new complex[nfft * Nb];

	for (i = 0; i < Nw; i++)
		dshift[i] = complex(cos(2 * pi * da * i), sin(2 * pi * da * i));
	sfm = complex(cos(2 * pi * R * da), sin(2 * pi * R * da));
	for (k = 0; k < KN; k++)
	{
		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(p);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < hnfft; i++)
		{
			Kft[i] = b[i];
			if (i > 0) Kft[nfft - i] = Kft[i].conj();
			mk = b[i].real;cp1 = b[i].imag;
			CPY2[i] = CPY2[i] + mk * mk + cp1 * cp1;
		}

		for (j = 1; j < Nb; j++) // calculate complete mesh across one deltaF
		{

			for (i = 0; i < Nw; i++)
				a[i] = a[i] * dshift[i]; // apply shifting
										 //fftw_execute(p); 
#ifdef USEFFTW
			fftw_execute(p);
#else
			ffthelper(a, b, Nw);
#endif // USEFFTW	
			inc = j * nfft;
			for (i = 0; i < nfft; i++)
				Kft[i + inc] = b[i];
		}

		for (i = 0; i < Nw; i++)
			dshift[i] = dshift[i] * sfm;//account for time of the next block

		for (i = 0; i < hnfft; i++)
		{
			sfa = Kft[i];
			mri = i;
			ci = 0;
			mk = -1.0;
			inc = 0;
			sptr = &SCoCmplx[i*a2];
			for (j = 0; j < a2; j++)
			{
				if (ci == Nb) {
					ci = 0;
					mri = mri - 1;
					sfa = sfa * complex(cos(2 * pi * (k * R) * Nb * da), sin(2 * pi * (k * R) * Nb * (mk * da)));
					if (mri == 0) {
						mk = 1.0F;
						sfa = sfa.conj();
					}
					if (mri == -1) inc = nfft;
					if ((mri < 0) && (opt == 0)) break; //prevent corr with negative spectrum
				}
				ri = inc + mri;
				kfk = complex(Kft[ri + ci * nfft].real, mk * Kft[ri + ci * nfft].imag);

				*sptr = *sptr + sfa * kfk;
				sptr++;
				ci++;
				if (Cohr)
					SCoX2[i*a2 + j] = SCoX2[i*a2 + j] + kfk.real * kfk.real + kfk.imag * kfk.imag;
				//SCoCmplx[i + j * hnfft] = SCoCmplx[i + j * hnfft] + kfk * sfa;

			}
		}
	}

	if (Cohr)
	{
		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
			{
				SCo[i + j * hnfft] = (float)(SCoCmplx[i*a2 + j].amp() / sqrt(CPY2[i] * SCoX2[i*a2 + j])); // pick the nearest bin
			}
		}
	}
	else
	{
		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
				SCo[i + j * hnfft] = (float)(SCoCmplx[i*a2 + j].amp() / KMU);
		}
	}

	fftw_destroy_plan(p);
	fftw_free(a);
	fftw_free(b);
	delete[] windat; delete[] SCoCmplx; delete[] SCoX2;
	delete[] dshift; delete[] CPY2; delete[] Kft;
	return 0;
}
int SCoh(float* x, float da, int a1, int a2, int L, int Noverlap, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using ACP

	da: delta alpha relative to Fs
	x() input data
	a1: minimum cyclic freq, must be 0, this is ignored input just to keep compatibility with other functions form
	a2: maximum cyclic frequency count such that max. alpha = da * a2
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data: of size Nw/2 * a2
	*/
	int nfft, hnfft, KN, i, k, j, R;
	R = Nw - Noverlap;
	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out;
	fftw_plan p;

	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	p = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw;i++)
		windat[i] = 2.0 * (0.5 - 0.5 * cos(2 * i * pi / Nw));

	KN = (int)((L - Noverlap) / (double)(Nw - Noverlap));
	//std::cout << "KN: " << KN << "\n";
	auto yws = new complex[hnfft * KN];
	auto dshift = new complex[L];
	//complex sfa, sfm;
	auto CPY2 = new double[hnfft];
	auto CPX2 = new double[hnfft];
	auto CPS = new complex[hnfft];
	auto xc = new complex[L];
	double KMU = 0.0;
	for (i = 0; i < L; i++) {
		xc[i] = complex(x[i], 0.0);
		dshift[i] = complex(cos(2 * pi * da * i), sin(2 * pi * da * i));
	}
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = KN * KMU;

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(p);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < hnfft; i++)
		{
			yws[i + k * hnfft] = b[i];
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
	}


	for (j = 0; j < a2; j++)
	{
		for (i = 0; i < hnfft; i++)
		{
			CPX2[i] = 0.0;
			CPS[i] = complex(0.0, 0.0);
		}
		for (k = 0; k < KN; k++) {
			for (i = 0; i < nfft; i++)
				a[i] = windat[i] * xc[i + k * R];
#ifdef USEFFTW
			fftw_execute(p);
#else
			ffthelper(a, b, Nw);
#endif // USEFFTW	
			for (i = 0; i < hnfft; i++)
			{
				if ((j * da <= ((float)i / Nw)) || (opt > 0)) {
					CPS[i] = CPS[i] + yws[i + k * hnfft] * (b[i].conj());
				}
				CPX2[i] = CPX2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
			}
		}
		for (i = 0; i < hnfft; i++) {
			if ((opt == 0) || (opt == 2)) {
				SCo[i + (j - a1) * hnfft] = (float)(CPS[i].amp() / KMU);
			}
			else {
				SCo[i + (j - a1) * hnfft] = (float)(CPS[i].amp() / sqrt((CPX2[i] * CPY2[i])));
			}

		}
		for (i = 0; i < L; i++)
			xc[i] = xc[i] * dshift[i];

	}



	fftw_destroy_plan(p);
	fftw_free(a);
	fftw_free(b);
	delete[] windat;
	delete[] yws; delete[] dshift; delete[] xc; delete[] CPY2; delete[] CPX2; delete[] CPS;
	return 0;
}

int AntoniFSC2(float* x, float &da, int a1, int a2, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using Antoni et al. Fast SC
	This is two sided, p = -P ... P (accprding to matlab code)
	da: delta alpha relative to Fs
	x() input data
	a1: minimum cyclic freq, must be 0, this is ignored input just to keep compatibility with other functions form
	a2: maximum cyclic frequency count such that max. alpha = da * a2
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, P, M, ir, ir2;
	float amax;
	complex sfi, sfj, sfj1;
	amax = a2 * da; //relative to Fs
	R = (int)floor(1.0f / (2.0f * amax));
	if (R > (Nw / 4))
		R = Nw / 4;
	//std::cout << "R: " << R << "\n";

	DbgMsg("AntoniFSC R: " << R << "\n");
	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out, *ink, *outk;
	fftw_plan fp, fp2;

	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw;i++)
		windat[i] = (0.5 - 0.5 * cos(2 * i * pi / Nw));
	P = (int)floor((float)Nw / (2.0f * R));
	//std::cout << "P: " << P << "\n";
	KN = (int)((L - Nw + R) / R);
	KNP2 = (int)pow(2.0, ceil(log2((float)KN)));//next power of two
												//std::cout << "KN, KNP2: " << KN<<", "<< KNP2 << "\n";
	DbgMsg("AntoniFSC KN, KNP2: " << KN << ", " << KNP2 << "\n");
	da = 1.0F / (KNP2 * R); //updated normalized da
							//std::cout << "updated da: " << da << "\n";
	complex* Ak = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	complex* Bk = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	ink = (fftw_complex*)&Ak[0]; //create alias
	outk = (fftw_complex*)&Bk[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink, outk, FFTW_FORWARD, FFTW_ESTIMATE);
	//std::cout << "KN: " << KN << "\n";
	auto yws = new complex[hnfft * KN];
	auto SCoCmplx = new complex[hnfft * a2];
	auto CPY2 = new double[hnfft];
	auto Rw = new double[KNP2]();
	auto SRw = new double[KN / 2]();
	double KMU = 0.0, RW0, w;

	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	RW0 = KMU;
	KMU = KN * KMU;
	M = Nw / 2 - 1;
	for (i = 0; i < KNP2;i++) //generate Rw(alpha) from 0 to KN-1
	{
		for (j = 0; j < Nw; j++)
			Rw[i] = Rw[i] + windat[j] * windat[j] * cos(2.0f * pi * (j - M) * i * da);
	}

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(fp);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < hnfft; i++)
		{
			yws[i + k * hnfft] = b[i];
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
	}
	if ((opt == 1) || (opt == 3)) {
		for (k = 0;k < KN;k++) {
			for (i = 0;i < hnfft;i++)
				yws[i + k * hnfft] = yws[i + k * hnfft] / (sqrt(CPY2[i] / KMU));
		}
	}
	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}
	M = KN / 2;
	if (a2 < M) M = a2;
	w = 2.0f * pi * (Nw / 2 - 1);
	sfi = complex(cos(w * da), -sin(w * da));
	for (j = 0; j <= P; j++) {
		sfj1 = complex(cos(w * j / Nw), sin(w * j / Nw));
		for (k = 0;k < hnfft;k++) {
			if (k >= j) {
				for (i = 0; i < KN; i++) {
					Ak[i] = yws[k + i * hnfft] * yws[(k - j) + i * hnfft].conj();
				}
#ifdef USEFFTW
				fftw_execute(fp2);
#else
				ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	

				sfj = sfj1;
				for (i = 0; i < M; i++) {
					//Bk[i] = Bk[i] * complex(cos(w * (i * da - (float)j / Nw)), -sin(w * (i * da - (float)j / Nw)));
					SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i] * sfj;
					sfj = sfj * sfi;
				}
				//now will process p < 0
				if (j > 0) {
					for (i = 0; i < KN; i++) {
						Ak[i] = yws[k + i * hnfft].conj() * yws[(k - j) + i * hnfft];
					}
#ifdef USEFFTW
					fftw_execute(fp2);
#else
					ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
					sfj = sfj1.conj();
					for (i = 0; i < M; i++) {
						//Bk[i] = Bk[i] * complex(cos(w * (i * da + (float)j / Nw)), -sin(w * (i * da + (float)j / Nw)));
						SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i] * sfj;
						sfj = sfj * sfi;
					}
				}

			}
			else {//aliased (for coherence)
				if (opt > 1) {
					for (i = 0; i < KN; i++) {
						Ak[i] = yws[k + i * hnfft] * yws[(j - k) + i * hnfft];
					}
#ifdef USEFFTW
					fftw_execute(fp2);
#else
					ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
					sfj = sfj1.conj();
					for (i = 0; i < M; i++) {
						//Bk[i] = Bk[i] * complex(cos(w * (i * da + (float)j / Nw)), -sin(w * (i * da + (float)j / Nw)));
						Bk[i] = Bk[i] * sfj;
						sfj = sfj * sfi;
						SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i];
					}
				}
			}
		}
	}
	for (i = 0; i < M; i++)
		SRw[i] = Rw[i]; //Rw(alpha - 0)
	for (j = 1; j <= P;j++) {
		for (i = 0; i < M; i++) {
			ir = i - (int)(j / (da * Nw)); //calculate shifted Rw(alpha - p * deltaF)
			ir2 = i + (int)(j / (da * Nw)); //calculate shifted Rw(alpha - p * deltaF) when p < 0
			if (ir < 0)  ir = -ir; //Rw(alpha) is even function
			SRw[i] = Rw[ir] + Rw[ir2] + SRw[i];
		}
	}
	for (i = 0; i < hnfft; i++)
	{
		for (j = 0; j < M; j++) {
			SCo[i + j * hnfft] = RW0 * SCoCmplx[i * a2 + j].amp() / KMU / SRw[j];
		}

	}
	fftw_destroy_plan(fp);
	fftw_free(a);
	fftw_free(b);
	fftw_destroy_plan(fp2);
	fftw_free(Ak);
	fftw_free(Bk);
	delete[] windat;
	delete[] yws; delete[] Rw; delete[] SRw; delete[] CPY2;
	delete[] SCoCmplx;
	return 0;
}


int FAM(float* x, float &da, int a1, int a2, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using using FFT Accumulation Method, Roberts et al.

	da: delta alpha relative to Fs
	x(): input data
	a1: minimum cyclic freq, must be 0, this is ignored input just to keep compatibility with other functions form
	a2: maximum cyclic frequency count such that max. alpha = da * a2
	L: input data length
	to get best speed, make sure L gives power of two for KN , where KN = (L - Nw + R) / R;=> L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data: size Nw*a2
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, ap, i1, j1, q, qr, R0, k1, aindex, findex, mf;
	float alpha0;

	fftw_complex *in, *out, *ink2, *outk2;
	complex *a, *b, *Ak, *Bk;
	fftw_plan fp, fp2;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nw); //new complex[Nw];
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nw); //new complex[Nw];
	a = (complex*)&in[0];
	b = (complex*)&out[0];

	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	auto windat = new double[Nw];


	nfft = Nw;
	hnfft = nfft / 2;
	//Prep WinData;
	for (i = 0; i < Nw;i++)
		windat[i] = (0.54 - 0.46 * cos(2 * i * pi / Nw)); //hamming window
	R = Nw / 4; //initial
	R0 = R;
	do {
		KN = (int)((L - Nw + R) / R);
		KNP2 = (int)pow(2, ceil(log2(KN)));//next power of two	
		DbgMsg("FAM R: " << R << ", KN: " << KN << ", KNP2: " << KNP2 << "\n";);
		R = (int)round((float)R * KN / KNP2);//try to reduce R to get best power-of-two KN		
		if (R == R0) break;
		R0 = R;
	} while (abs(KNP2 - KN) > 20);
	qr = KNP2 *  R / (2 * Nw);

	ink2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * KNP2); //new complex[KNP2];
	outk2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * KNP2); //new complex[KNP2];
	Ak = (complex*)&ink2[0]; //create complex alias
	Bk = (complex*)&outk2[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink2, outk2, FFTW_FORWARD, FFTW_ESTIMATE);

	da = 1.0F / (KNP2 * R); //updated normalized da
							//std::cout << "KN, KNP2: " << KN << ", " << KNP2 << "\n";
							//std::cout << "updated da: " << da << "\n";
	auto yws = new complex[nfft * KN];
	auto SCoCmplx = new complex[(nfft + 1) * a2];
	auto CPY2 = new double[nfft];
	double KMU = 0.0, w;

	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];
	KMU = KN * KMU / 2.0f;


	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	w = 2.0f * pi * R / Nw;
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(fp);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	

		for (i = 0; i < nfft; i++)
		{
			yws[i + k * nfft] = b[i] * complex(cos(w * i * k), -sin(w * i * k)); //phase compensation;
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
	}
	if ((opt == 1) || (opt == 3)) {
		for (k = 0;k < KN;k++) {
			for (i = 0;i < nfft;i++)
				yws[i + k * nfft] = yws[i + k * nfft] / (sqrt(CPY2[i] / KMU));//to calculate coherence
		}
	}
	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}

	mf = nfft;
	if (opt < 1) mf = hnfft;
	for (j = 0; j < mf; j++) {
		for (i = 0;i < mf; i++) {

			for (k = 0; k < KN; k++) {
				Ak[k] = yws[j + k * nfft] * yws[i + k * nfft].conj();
			}
#ifdef USEFFTW
			fftw_execute(fp2);
#else
			ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
			if (i >= hnfft) {
				i1 = i - nfft;
			}
			else {
				i1 = i;
			}
			if (j >= hnfft) {
				j1 = j - nfft;
			}
			else {
				j1 = j;
			}
			//f0 = (i1 + j1) / 2;
			findex = i1 + j1;
			alpha0 = (float)(j1 - i1) / (Nw * da);
			//
			ap = (int)floor(alpha0);
			//DbgMsg(" ap: " << ap);
			for (q = -qr; q <= qr; q++) {
				aindex = ap + q;
				if ((aindex >= 0) && (findex >= 0) && (aindex < a2)) {
					if (q < 0) {
						k1 = KNP2 + q;
					}
					else {
						k1 = q;
					}

					SCoCmplx[findex * a2 + aindex] = Bk[k1] + SCoCmplx[findex * a2 + aindex];
				}

			}

		}
	}

	for (i = 0; i < nfft; i++)
	{
		for (j = 0; j < a2; j++) {
			SCo[i + j * nfft] = (float)(SCoCmplx[i * a2 + j].amp() / KMU);
		}

	}

	fftw_destroy_plan(fp);
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(fp2);
	fftw_free(ink2);
	fftw_free(outk2);

	delete[] windat;
	delete[] yws; delete[] CPY2;
	delete[] SCoCmplx;
	return 0;
}

int FDirSC(float* x, float &da, int a1, int a2, int L, int Nw, int opt, float* SCo)
{
	/*
	calculate spectral correlation/coherance using  Faster SC, Borghesani and Antoni 2018

	da: delta alpha relative to Fs
	x() input data
	a1: minimum cyclic freq, must be 0, this is ignored input just to keep compatibility with other functions form
	a2: maximum cyclic frequency count such that max. alpha = da * a2
	L: input data length
	to get best speed and accuracy, make sure L gives power of two for KN , where KN = (L - Nw + R) / R; => L = R*KN + Nw - R
	Nw: window length        '
	opt:0 :SC no negative frequencies correlation, 1 : Scoh, 2: SC with negative frequencies, 3: Scoh with negative freq

	SCo() output data
	*/
	int nfft, hnfft, KN, KNP2, i, k, j, R, P, M, ri, KNR;
	float amax, Nb;
	amax = a2 * da; //relative to Fs
	R = (int)floor(1.0f / (2.0f * amax));
	//std::cout << "FDirSC, R: " << R << "\n";
	if (R > (Nw / 4)) R = Nw / 4;
	DbgMsg("FDirSC, R: " << R << "\n");

	nfft = Nw;
	hnfft = nfft / 2;
	auto windat = new double[Nw];
	fftw_complex *in, *out, *ink, *outk, *inr;
	fftw_plan fp, fp2, fp3;

	complex* a = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	complex* b = (complex*)fftw_malloc(sizeof(complex) * Nw); //new complex[Nw];
	in = (fftw_complex*)&a[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)&b[0];
	fp = fftw_plan_dft_1d(Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//Prep WinData;
	for (i = 0; i < Nw;i++)
		windat[i] = (0.5 - 0.5 * cos(2 * i * pi / Nw));
	P = (int)floor((float)Nw / (2.0f * R));
	//std::cout << "P: " << P << "\n";
	KN = (int)((L - Nw + R) / R);
	KNP2 = (int)pow(2.0, ceil(log2(KN)));//next power of two
										 //std::cout << "KN, KNP2: " << KN << ", " << KNP2 << "\n";
	DbgMsg("FDirSC KN, KNP2: " << KN << ", " << KNP2 << "\n");
	da = 1.0F / (KNP2 * R); //updated normalized da
							//std::cout << "updated da: " << da << "\n";
	Nb = (1.0f / da) / Nw;
	complex* Ak = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	complex* Bk = (complex*)fftw_malloc(sizeof(complex) * KNP2); //new complex[Nw];
	ink = (fftw_complex*)&Ak[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	outk = (fftw_complex*)&Bk[0];
	fp2 = fftw_plan_dft_1d(KNP2, ink, outk, FFTW_FORWARD, FFTW_ESTIMATE);
	//std::cout << "KN: " << KN << "\n";
	auto yws = new complex[hnfft * KN];
	auto ywsd = new complex[hnfft * KN];
	auto SCoCmplx = new complex[hnfft * a2];
	auto CPY2 = new double[hnfft];
	auto Dp = new complex[Nw];
	double KMU = 0.0f, w;
	complex s1;
	KNR = KNP2*R;
	complex* Ar = (complex*)fftw_malloc(sizeof(complex) * KNR); //new complex[Nw];	
	inr = (fftw_complex*)&Ar[0]; //(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	fp3 = fftw_plan_dft_1d(KNR, inr, inr, FFTW_FORWARD, FFTW_ESTIMATE);//in-place transform
	for (i = 0; i < Nw; i++)
		KMU = KMU + windat[i] * windat[i];

	M = Nw / 2 - 1;
	w = 2.0f * pi / Nw;
	for (i = 0; i < Nw;i++) //generate Dirichlet kernel
	{
		s1 = complex(0, 0);
		for (j = -P; j <= P; j++)
			s1 = s1 + complex(cos(w * j * (i - hnfft)), sin(w * j * (i - hnfft)));
		Dp[i] = s1;
	}
	for (i = 0; i < KNR;i++) //generate calibration data
	{
		if (i < Nw) {
			Ar[i] = windat[i] * windat[i] * Dp[i];
		}
		else {
			Ar[i] = complex(0.0f, 0.0f);
		}
	}
	fftw_execute(fp3); //calibration data now in Ar[]

	for (i = 0; i < hnfft; i++)
	{
		CPY2[i] = 0.0;
	}
	for (k = 0; k < KN; k++) // calculate initial spectral data
	{

		for (i = 0; i < nfft; i++)
			a[i] = complex(x[i + k * R] * windat[i], 0.0);
#ifdef USEFFTW
		fftw_execute(fp);
#else
		ffthelper(a, b, Nw);
#endif // USEFFTW	
		for (i = 0; i < hnfft; i++)
		{
			yws[i + k * hnfft] = b[i];
			CPY2[i] = CPY2[i] + b[i].real * b[i].real + b[i].imag * b[i].imag;
		}
		for (i = 0; i < nfft; i++)
			a[i] = a[i] * Dp[i];
		//fftw_execute(fp);
		ffthelper(a, b, Nw);
		for (i = 0; i < hnfft; i++)
		{
			ywsd[i + k * hnfft] = b[i];
		}
	}

	for (k = 0; k < KNP2; k++) {
		Ak[k] = complex(0.0f, 0.0f);
	}
	M = KN / 2;
	if (a2 < M) M = a2;
	for (k = 0; k < hnfft; k++) {
		for (i = 0; i < KN; i++) {
			Ak[i] = yws[k + i * hnfft] * ywsd[k + i * hnfft].conj();
		}
#ifdef USEFFTW
		fftw_execute(fp2);
#else
		ffthelper(Ak, Bk, KNP2);
#endif // USEFFTW	
		for (i = 0; i < M; i++)
			SCoCmplx[k * a2 + i] = SCoCmplx[k * a2 + i] + Bk[i] / Ar[i];

	}

	if ((opt == 0) || (opt == 2))
	{

		for (i = 0; i < hnfft; i++)
		{
			for (j = 0; j < a2; j++)
				SCo[i + j * hnfft] = (float)(SCoCmplx[i*a2 + j].amp() / KN);
		}
	}
	else
	{
		for (i = 0; i < hnfft; i++)
		{

			for (j = 0; j < a2; j++)
			{
				ri = (int)round(i - (j / Nb));
				if (ri < 0)
				{
					ri = -ri;
				}
				SCo[i + j * hnfft] = (float)(KMU * SCoCmplx[i*a2 + j].amp() / sqrt(CPY2[i] * CPY2[ri])); // pick the nearest bin
			}
		}
	}
	fftw_destroy_plan(fp);
	fftw_free(a);
	fftw_free(b);
	fftw_destroy_plan(fp2);
	fftw_free(Ak);
	fftw_free(Bk);
	fftw_destroy_plan(fp3);
	fftw_free(Ar);

	delete[] windat;
	delete[] yws; delete[] ywsd; delete[] Dp; delete[] CPY2;
	delete[] SCoCmplx;
	return 0;
}



int main()
{
	float da, wdf,SamRate,amax;
	int a1, a2,L, Nw, Noverlap, i,ia;
	SamRate = 32768;
	
	srand(time(NULL));   // Initialization, should only be called once.
	
	L = 131072;
	
	Nw = 512;	
	Noverlap = (int)(Nw * 0.6667);//or R = 0.333Nw for ACP and Fast ACP
	
	wdf = SamRate / Nw;//delta freq
	da = (wdf / 64) / SamRate; //relative to Sample rate
	auto SCo = new float[Nw  * L/2];
	auto x = new float[L+32]();//putting () will initilize to zero all elements
	a1 = 0;
	a2 = 4096;
	
	amax = 0.03125;
	for (ia = 1;ia < 6;ia++) {
		a2 = (int)(amax*L);
		std::cout << "amax, alpha count: " << amax <<", " << a2 << "\n";
		

		for (i = 0;i < L;i++) {
			x[i] = (1.0f+sin(2.0 * pi * i * 10 / Nw)) * (float)rand() / 32767.0f;
		}

		auto start = std::chrono::high_resolution_clock::now();
		FSCoh(x, da, a1, a2, L, Noverlap, Nw, 0, SCo);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = finish - start;
		std::cout << "Elapsed time fast ACP: " << elapsed1.count() << " s\n";

		SCoh(x, da, a1, a2, L, Noverlap, Nw, 0, SCo);
		auto finish2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed2 = finish2 - finish;
		std::cout << "Elapsed time normal ACP: " << elapsed2.count() << " s\n";

		AntoniFSC2(x, da, a1, a2, L, Nw, 0, SCo);
		auto finish3 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed3 = finish3 - finish2;
		std::cout << "Elapsed time Antoni FSC2: " << elapsed3.count() << " s\n";

		FAM(x, da, a1, a2, L, Nw, 2, SCo);
		auto finish4 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed4 = finish4 - finish3;
		std::cout << "Elapsed time FAM: " << elapsed4.count() << " s\n";

		FDirSC(x, da, a1, a2, L, Nw, 0, SCo);
		auto finish5 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed5 = finish5 - finish4;
		std::cout << "Elapsed time Fast Dir SC: " << elapsed5.count() << " s\n";

		amax = 2.0f*amax;
		if (amax > 0.5) break;
		std::cout << "============================================================ \n";
	}
	delete[] x; delete[] SCo;
	getchar();
    return 0;
}

