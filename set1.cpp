#include "Source.h"
#include "parameter.h"

using namespace std;

//default_random_engine e(199161);
//
//uniform_real_distribution<> u(-1, 1);


int main()
{
	clock_t init, final;
	double MTeach = 7;
	init = clock();
	int i, j, k, t, q,i1,j1, i2, j2, k2, i3, j3, k3;

	//	const int n = 200,m=20;
	double Dc = Dc0;// 1-10 typical Cell Signalling : Changing Shape Changes the Signal
	double vm = vm0;// Tracking Single Kinesin Molecules in the Cytoplasm of Mammalian Cells
	double km = km0;
	double kf = kf0;//Tracking Single Kinesin Molecules in the Cytoplasm of Mammalian Cells
	double ka = ka0;//
	double ka1 = ka10;//
	double kd = kd0;// 0.1-100.... Cell Signalling : Changing Shape Changes the Signal
	double kp = kp0;//review paper
	double kp1 = kp10;
	double KMT = KMT0;
	double k0 = k00;//old paper
	double kc = kc0;//review
	double Kca = Kca0;



	char output_L[120];
	char output_temp[120];
	sprintf(output_L, "L%d/", L);
	/*_mkdir(output_L);*/
	FILE *time;
	sprintf(output_temp, "%stime.dat", output_L);
	time = fopen(output_temp, "w");
	FILE *kemoCmu;
	sprintf(output_temp, "%sCmua.dat", output_L);
	kemoCmu = fopen(output_temp, "w");
	FILE *kemoCa;
	sprintf(output_temp, "%sCaa.dat", output_L);
	kemoCa = fopen(output_temp, "w");
	FILE *kemoCmd;
	sprintf(output_temp, "%sCmda.dat", output_L);
	kemoCmd = fopen(output_temp, "w");
	FILE* outCa;
	sprintf(output_temp, "%sCa.dat", output_L);
	outCa = fopen(output_temp, "w");
	FILE* outCc;
	sprintf(output_temp, "%sCc.dat", output_L);
	outCc = fopen(output_temp, "w");
	FILE* out;
	sprintf(output_temp, "%sCcsum.dat", output_L);
	out = fopen(output_temp, "w");
	FILE* outMTu;
	sprintf(output_temp, "%sMTu.dat", output_L);
	outMTu = fopen(output_temp, "w");
	FILE* outMTd;
	sprintf(output_temp, "%sMTd.dat", output_L);
	outMTd = fopen(output_temp, "w");
	FILE* outMT0;
	sprintf(output_temp, "%sMT0.dat", output_L);
	outMT0 = fopen(output_temp, "w");
	FILE* outMT;
	sprintf(output_temp, "%sMT.dat", output_L);
	outMT = fopen(output_temp, "w");
	FILE* outCmd;
	sprintf(output_temp, "%sCmd.dat", output_L);
	outCmd = fopen(output_temp, "w");
	FILE* outCmu;
	sprintf(output_temp, "%sCmu.dat", output_L);
	outCmu = fopen(output_temp, "w");
	FILE* outCau;
	sprintf(output_temp, "%sCau.dat", output_L);
	outCau = fopen(output_temp, "w");
	FILE* outCad;
	sprintf(output_temp, "%sCad.dat", output_L);
	outCad = fopen(output_temp, "w");
	FILE* outCar;
	sprintf(output_temp, "%sCar.dat", output_L);
	outCar = fopen(output_temp, "w");
	FILE* xdir;
	sprintf(output_temp, "%sxdir.dat", output_L);
	xdir = fopen(output_temp, "w");
	FILE* testratesca;
	sprintf(output_temp, "%stestrates.dat", output_L);
	testratesca = fopen(output_temp, "w");
	FILE* testratesMT;
	sprintf(output_temp, "%stestratesMT.dat", output_L);
	testratesMT = fopen(output_temp, "w");
	FILE* kemoCc;
	sprintf(output_temp, "%sCckymo.dat", output_L);
	kemoCc = fopen(output_temp, "w");


	FILE *para;
	sprintf(output_temp, "%sparameter.dat", output_L);
	para = fopen(output_temp, "w");


	FILE *pa;
	sprintf(output_temp, "%srates.dat", output_L);
	pa = fopen(output_temp, "w");


	/*snprintf(output_L, 120,".\\%d\\", L);
	mkdir(output_L, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);



	FILE *time;
	snprintf(output_temp,120, "%stime.dat", output_L);
	time = fopen(output_temp, "w");
	FILE *kemoCmu;
	snprintf(output_temp, 120, "%sCmua.dat", output_L);
	kemoCmu = fopen(output_temp, "w");
	FILE *kemoCa;
	snprintf(output_temp, 120, "%sCaa.dat", output_L);
	kemoCa = fopen(output_temp, "w");
	FILE *kemoCmd;
	snprintf(output_temp, 120, "%sCmda.dat", output_L);
	kemoCmd = fopen(output_temp, "w");
	FILE* outCa;
	snprintf(output_temp, 120, "%sCa.dat", output_L);
	outCa = fopen(output_temp, "w");
	FILE* outCc;
	snprintf(output_temp, 120, "%sCc.dat", output_L);
	outCc = fopen(output_temp, "w");
	FILE* out;
	snprintf(output_temp, 120, "%sCcsum.dat", output_L);
	out = fopen(output_temp, "w");
	FILE* outMTu;
	snprintf(output_temp, 120, "%sMTu.dat", output_L);
	outMTu = fopen(output_temp, "w");
	FILE* outMTd;
	snprintf(output_temp, 120, "%sMTd.dat", output_L);
	outMTd = fopen(output_temp, "w");
	FILE* outMT0;
	snprintf(output_temp, 120, "%sMT0.dat", output_L);
	outMT0 = fopen(output_temp, "w");
	FILE* outMT;
	snprintf(output_temp, 120, "%sMT.dat", output_L);
	outMT = fopen(output_temp, "w");
	FILE* outCmd;
	snprintf(output_temp, 120, "%sCmd.dat", output_L);
	outCmd = fopen(output_temp, "w");
	FILE* outCmu;
	snprintf(output_temp, 120, "%sCmu.dat", output_L);
	outCmu = fopen(output_temp, "w");
	FILE* outCau;
	snprintf(output_temp, 120, "%sCau.dat", output_L);
	outCau = fopen(output_temp, "w");
	FILE* outCad;
	snprintf(output_temp, 120, "%sCad.dat", output_L);
	outCad = fopen(output_temp, "w");
	FILE* outCar;
	snprintf(output_temp, 120, "%sCar.dat", output_L);
	outCar = fopen(output_temp, "w");
	FILE* xdir;
	snprintf(output_temp, 120, "%sxdir.dat", output_L);
	xdir = fopen(output_temp, "w");
	FILE* testratesca;
	snprintf(output_temp, 120, "%stestrates.dat", output_L);
	testratesca = fopen(output_temp, "w");
	FILE* testratesMT;
	snprintf(output_temp, 120, "%stestratesMT.dat", output_L);
	testratesMT = fopen(output_temp, "w");
	FILE* kemoCc;
	snprintf(output_temp, 120, "%sCckymo.dat", output_L);
	kemoCc = fopen(output_temp, "w");*/

	double Ca[ai], Catemp[ai];
	double Csum;
	double(*Cc)[r + 1], (*Cctemp)[r + 1];
	double(*Cmu)[r + 1], (*Cmutemp)[r + 1], (*Cmd)[r + 1], (*Cmdtemp)[r + 1], (*F)[r + 1];
	double dCc, dCmu, dCmd;//differentation of variables
	double MTu[r + 1], MTd[r + 1], MT0[r + 1];//global variables
	double  dMTu[r + 1], dMTd[r + 1], dMT0[r + 1];//differentation of global variables
	double MTutemp[r + 1], MTdtemp[r + 1], MT0temp[r + 1];//temp global variables
	double Kau[r + 1], Kad[r + 1], Kmtu[r + 1], Kmtd[r + 1]; //hill equation values
	double A[x + 1], B[x + 1], C[r + 1], D[r + 1], Fx[x + 1], Fr[x + 1], rx[x + 1], rr[r + 1], Ccx[x + 1], Ccr[r + 1], Hx[x + 1], Hr[r + 1], A1[x + 1], C1[r + 1];
	double ar[r], dr[r];
	double  dx[x];
	double Ccsum[x + 1];
	double diff;
	double conu;
	double cond;
	double Casum = 0;
	double MTsum = 0, Cmusum = 0, Cmdsum = 0;
	double MTus = 0;
	double MTds = 0;
	double MT0s = 0, Ccsum1 = 0;
	double Cmds = 0;
	double Cmus = 0,Ccs=0;
	double Cavau[x + 1], Cavad[x + 1];
	Cc = new double[x + 1][r + 1];
	Cctemp = new double[x + 1][r + 1];
	Cmu = new double[x + 1][r + 1];
	Cmutemp = new double[x + 1][r + 1];
	Cmd = new double[x + 1][r + 1];
	Cmdtemp = new double[x + 1][r + 1];
	F = new double[x + 1][r + 1];
	double rn[r + 1];
	formdx(dx);
	formdr(dr);
	rn[0] = 0;
	for (j = 1; j < r + 1; j++)
	{
		rn[j] = rn[j - 1] + dr[j - 1];
		//cout << j << "\t" << rn[j] << "\t";
	}

	double R = 4 * rad / ((2 * rad - dr[r - 1])*dr[r - 1]);
	double X = 2 / dx[0];

	double xd = 0;
	
	for (i = 0; i < x; i++)
	{
			fprintf(xdir, "%lf\t", xd);
			xd = xd + dx[i];

	}
	fclose(xdir);
	double Cau = 0, Cad = 0;
	t = tf / dt + 1;

	//This section is initialze 
	for (q = 0; q < ai; q++)
	{
		Ca[q] = 0;
		Catemp[q] = 0;
	}
	for (i = 0; i < x + 1; i++)
	{
		for (j = 0; j < r + 1; j++)
		{
			F[i][j] = 0;
			Cc[i][j] = 0;
			Cctemp[i][j] = 0;
			Cmu[i][j] = 0;
			Cmutemp[i][j] = 0;
			Cmd[i][j] = 0;
			Cmdtemp[i][j] = 0;
		}
	}
	for (j = 0; j < r + 1; j++)
	{
		Hr[j] = 0;
		Fr[j] = 0;
		rr[j] = 0;
		Ccr[j] = 0;
		Kau[j] = 0;
		Kad[j] = 0;
		MTu[j] = 0;
		MTd[j] = 0;
		MT0[j] = 0;
		MTutemp[j] = 0;
		MTdtemp[j] = 0;
		MT0temp[j] = 0;
		dMTu[j] = 0;
		dMTd[j] = 0;
		dMT0[j] = 0;
	}

	for (i = 0; i < x + 1; i++)
	{
		Cavau[i] = 0;
		Cavad[i] = 0;
		Ccsum[i] = 0;
		Hx[i] = 0;
		Fx[i] = 0;
		rx[i] = 0;
		Ccx[i] = 0;
	}
	//done initialize


	//FILE * input;
	//input = fopen("Ca1.dat", "r");
	//for (q = 0; q < ai; q++)
	//{
	//	fscanf(input, "%lf\t", &Ca[q]);
	//}
	double V = 0;
	//for (i = 1; i <x; i++)
	//{
	//	Ca[2 * r+i] = 100 / (2*pi*rad*L);
	//}

	//for (i = 0; i <x + 1; i++)
	//{
	//	Ca[2 * r + i] = 100 / (2 * pi*rad*L);
	//}
	for (i =23; i < 29; i++)
	{
		for (j = 0; j < r; j++)
		{
			V = V + dx[i] * 2 * pi*(rn[j] + dr[j] / 2)*dr[j];
		}
	}

	double Cceach = Ccall/ V;
	for (i = 23; i <29; i++)
	{
		for (j = 0; j < r + 1; j++)
		{
			Cc[i][j] = Cceach;
		}
	}
	//for (j = r - 3; j < r; j++)
	//{
	//	for (i = 25; i < 30; i++)
	//	{
	//		V = V + dx[i] * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] / 2 + dx[i] * 2 * pi*(rn[j - 1] + dr[j - 1] / 2)*dr[j - 1] / 2;
	//	}
	//}
	//double Cceach = 100 / V;
	//for (j = r - 3; j < r; j++)
	//{
	//	for (i = 25; i < 30; i++)
	//	{
	//		Cc[i][j] = Cceach;
	//	}
	//}

	/*for (j = r - 5; j < r; j++)
	{
		for (i = 25; i < 30; i++)
		{
			V = V + dx[i] * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] / 2 + dx[i] * 2 * pi*(rn[j - 1] + dr[j - 1] / 2)*dr[j - 1] / 2;
		}
	}
	double Cceach = 100 / V;
	for (j = r - 5; j < r; j++)
	{
		for (i = 25; i < 30; i++)
		{
			Cc[i][j] = Cceach;
		}
	}*/

	for (j = 0; j < r + 1; j++)
	{
		MTu[j] = MTeach;
		MTutemp[j] =MTeach;
		MTd[j] = 1;
		MTdtemp[j] =1;
		MT0[j] = 1;
		MT0temp[j] = 1;
	}



	Csum = 0;
	k = 0;
	//char buffer[50];



	//sprintf(buffer, "Z%d.dat", k);
	//ofstream xxfile(buffer, ios::out);

	for (i = 0; i < x; i++)
	{

		for (j = 0; j < r; j++)
		{
			Ccsum[i] = Ccsum[i] + 0.25*(Cc[i][j] + Cc[i + 1][j] + Cc[i][j + 1] + Cc[i + 1][j + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] * dx[i];

		}
		//xxfile << Ccsum[i] << "\t";
		Csum = Csum + Ccsum[i];
	}

	//xxfile << Ccsum[x] << "\t";
	Csum = Csum + Ccsum[x];
	//xxfile.close();
	//cout << "Sum from Cc\t" << Csum << endl;
	//cin.get();

	Casum = 0;
	for (j = 0; j < r - 1; j++)
	{
		Casum = Casum + 0.5*(Ca[2 * j] + Ca[2 * (j + 1)] + Ca[2 * j + 1] + Ca[2 * (j + 1) + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j];
	}
	Casum = Casum + 0.5*(Ca[2 * (r - 1)] + Ca[2 * r] + Ca[2 * r - 1] + Ca[2 * r + x]) * 2 * pi*(rn[j] + dr[r - 1] / 2)*dr[r - 1];
	//x line
	for (i = 0; i < x; i++)
	{
		Casum = Casum + 0.5*(Ca[2 * r + i] + Ca[2 * r + i + 1]) * 2 * pi*rad*dx[i];
	}

	Csum = Csum + Casum + Cmusum + Cmdsum;

	//fprintf(out, "\n");
	cout << "Sum from Ca\t" << Casum << endl;

	formABA1(A, B, A1, dx, Dc);

	formCDC1(C, D, C1, dr, Dc);
	//from inside the k loop no noise simulation
	int kkk = 0;
	//......................................................
	for (k = 0; k < t; k++)
	{

		//Dc = Dc0*(1 + 0.01*u(e));
		//vm = vm0*(1 + 0.01*u(e));
		//km = km0*(1 + 0.01*u(e));
		//kf = kf0*(1 + 0.01*u(e));
		//ka = ka0*(1 + 0.01*u(e));
		//ka1 = ka10*(1 + 0.01*u(e));
		//kd = kd0*(1 + 0.01*u(e));
		//kp = kp0*(1 + 0.01*u(e));
		//kp1 = kp10*(1 + 0.01*u(e));
		//KMT = KMT0*(1 + 0.01*u(e));
		//k0 = k00*(1 + 0.01*u(e));
		//kc = kc0*(1 + 0.01*u(e));//1% noise
		//Kca = Kca0*(1 + 0.01*u(e));

		for (q = 0; q < ai; q++)
		{
			Catemp[q] = Ca[q];
		}
		Csum = 0;
		for (i = 0; i < x + 1; i++)
		{
			for (j = 0; j < r + 1; j++)
			{
				Cctemp[i][j] = Cc[i][j];
				Cmutemp[i][j] = Cmu[i][j];
				Cmdtemp[i][j] = Cmd[i][j];
			}
		}
		for (j = 0; j < r + 1; j++)
		{
			MTutemp[j] = MTu[j];
			MTdtemp[j] = MTd[j];
			MT0temp[j] = MT0[j];
		}
		//make a copy

		//find F first

		for (j = 0; j < r; j++)
		{
			Kau[j] = FKa(Catemp[2 * j + 1],kc,Kca);
			Kad[j] = FKa(Catemp[2 * j], kc, Kca);
		}
		Kau[r] = FKa(Catemp[2 * r + x], kc, Kca);
		Kad[r] = FKa(Catemp[2 * r], kc, Kca);

		for (i = 0; i < x + 1; i++)
		{
			for (j = 0; j < r + 1; j++)
			{
				F[i][j] = -km*((MTutemp[j] - Cmutemp[i][j]) + (MTdtemp[j] - Cmdtemp[i][j]))*Cctemp[i][j] + kf*(Cmutemp[i][j] + Cmdtemp[i][j])
					+ Kau[j] * Cmutemp[i][j] + Kad[j] * Cmdtemp[i][j] + Cmutemp[i][j] / MTutemp[j] * k0*(MTutemp[j] - Cmutemp[x ][j])
					+ Cmdtemp[i][j] / MTdtemp[j] * k0*(MTdtemp[j] - Cmdtemp[0][j]);
			}
		}

		for (j = 0; j < r; j++)
		{
			Ca[2 * j] = Catemp[2 * j] + dt*FdCa(Cctemp[0][j], Catemp[2 * j], ka, ka1, kd);
			Ca[2 * j + 1] = Catemp[2 * j + 1] + dt*FdCa(Cctemp[x][j], Catemp[2 * j + 1], ka, ka1, kd);
		}
		for (i = 0; i < x + 1; i++)
		{
			Ca[2 * r + i] = Catemp[2 * r + i] + dt*FdCa(Cctemp[i][r], Catemp[2 * r + i], ka, ka1, kd);
		}//checked feb 28

		if (k % 2 == 0)//k is the time point
		{//x implicit r explicit

			for (i1 = 0; i1 < x + 1; i1++)//r=0
			{
				Hx[i1] = -(D[0] - 1 / dt)* Cctemp[i1][0] - C1[0] * Cctemp[i1][1]; 
				rx[i1] = Cctemp[i1][0] / dt + Hx[i1] + F[i1][0];//the other side of the equation
				Ccx[i1] = Cctemp[i1][0];//make a copy of Cc
			}

			rx[0] = Cctemp[0][0] / dt + Hx[0] - FdCa(Cctemp[0][0], Catemp[0], ka, ka1, kd)*X + F[0][0];//minus sign is derived from having a fake point that represented by flux and the point before
			rx[x] = Cctemp[x][0] / dt + Hx[x] - FdCa(Cctemp[x][0], Catemp[1], ka, ka1, kd)*X + F[x][0];//checked Mar 1

			tridag(A, B, A1, rx, Ccx, x + 1);

			for (i1 = 0; i1 < x + 1; i1++)
			{
				Cc[i1][0] = Ccx[i1];
			}

			for (j2 = 1; j2 < r; j2++)
			{//r from 1 to r-1
				for (i2 = 0; i2 < x + 1; i2++)
				{
					Hx[i2] = -Cctemp[i2][j2 - 1] * C[j2] - Cctemp[i2][j2] * (D[j2] - 1 / dt) - Cctemp[i2][j2 + 1] * C1[j2];
					rx[i2] = Cctemp[i2][j2] / dt + Hx[i2] + F[i2][j2];
					Ccx[i2] = Cctemp[i2][j2];
				}
				rx[0] = Cctemp[0][j2] / dt + Hx[0] - FdCa(Cctemp[0][j2], Catemp[2 * j2], ka, ka1, kd)*X + F[0][j2];//FdCa(Cctemp[0][j2], Catemp[2*j2]);
				rx[x] = Cctemp[x][j2] / dt + Hx[x] - FdCa(Cctemp[x][j2], Catemp[2 * j2 + 1], ka, ka1, kd) *X + F[x][j2];// FdCa(Cctemp[x - 1][j2], Catemp[2*j2+1]);
				tridag(A, B, A1, rx, Ccx, x + 1);
				for (i2 = 0; i2 < x + 1; i2++)
				{
					Cc[i2][j2] = Ccx[i2];
				}
			}
			//j=r
			for (i3 = 0; i3 < x + 1; i3++)
			{
				Hx[i3] = -C[r] * Cctemp[i3][r - 1] - (D[r] - 1 / dt)*Cctemp[i3][r];
				rx[i3] = Cctemp[i3][r] / dt + Hx[i3] - FdCa(Cctemp[i3][r], Catemp[2 * r + i3], ka, ka1, kd)*R + F[i3][r];
				Ccx[i3] = Cctemp[i3][r - 1];
			}
			rx[0] = Cctemp[0][r] / dt + Hx[0] - FdCa(Cctemp[0][r], Catemp[2 * r], ka, ka1, kd)*(R + X) + F[0][r];
			rx[x] = Cctemp[x][r] / dt + Hx[x] - FdCa(Cctemp[x][r], Catemp[2 * r + x], ka, ka1, kd)*(R + X) + F[x][r];

			tridag(A, B, A1, rx, Ccx, x + 1);
			for (i3 = 0; i3 < x + 1; i3++)
			{
				Cc[i3][r] = Ccx[i3];
			}

		}

		//============================

		else
		{//r implicit x explicit

			//x=0
			for (j1 = 0; j1 < r + 1; j1++)
			{
				Hr[j1] = -(B[0] - 1 / dt)*Cctemp[0][j1] - A1[0] * Cctemp[1][j1];
				rr[j1] = Cctemp[0][j1] / dt + Hr[j1] - FdCa(Cctemp[0][j1], Catemp[2 * j1], ka, ka1, kd)*X + F[0][j1];//+ F[ii][j]
				Ccr[j1] = Cctemp[0][j1];
			}
			rr[r] = Cctemp[0][r] / dt + Hr[r] - FdCa(Cctemp[0][r], Catemp[2 * r], ka, ka1, kd) *(R + X) + F[0][r];//FdCa(Cctemp[i2][r-1], Catemp[2 * (r-1) + i2]);
			tridag(C, D, C1, rr, Ccr, r + 1);
			for (j1 = 0; j1 < r + 1; j1++)
			{
				Cc[0][j1] = Ccr[j1];
			}

			//x=1->x
			for (i2 = 1; i2 < x; i2++)
			{
				for (j2 = 0; j2 < r + 1; j2++)
				{
					Hr[j2] = -Cctemp[i2 - 1][j2] * A[i2] - Cctemp[i2][j2] * (B[i2] - 1 / dt) - Cctemp[i2 + 1][j2] * A1[i2];
					rr[j2] = Cctemp[i2][j2] / dt + Hr[j2] + F[i2][j2];
					Ccr[j2] = Cctemp[i2][j2];
				}
				rr[r] = Cctemp[i2][r] / dt + Hr[r] - FdCa(Cctemp[i2][r], Catemp[2 * r + i2],ka,ka1,kd) *R+ F[i2][r];

				tridag(C, D, C1, rr, Ccr, r + 1);
				for (j2 = 0; j2 < r + 1; j2++)
				{
					Cc[i2][j2] = Ccr[j2];
				}
			}
			//i=x
			for (j3 = 0; j3 < r + 1; j3++)
			{
				Hr[j3] = -A[x] * Cctemp[x - 1][j3] - (B[x] - 1 / dt)* Cctemp[x][j3];//(Cctemp[i2 + 1][j] + Cctemp[i2 - 1][j] - 2 * Cctemp[i2][j]);//+ax*(Cctemp[ii + 1][j] - Cctemp[ii][j]) / (2 * j);;
				rr[j3] = Cctemp[x][j3] / dt + Hr[j3] - FdCa(Cctemp[x][j3], Catemp[2 * j3 + 1], ka, ka1, kd)*X + F[x][j3];//+ F[ii][j]
				Ccr[j3] = Cctemp[x][j3];
			}

			rr[r] = Cctemp[x][r] / dt + Hr[r] - FdCa(Cctemp[x][r], Catemp[2 * r + x], ka, ka1, kd)*(R + X) + F[x][r];//FdCa(Cctemp[i2][r-1], Catemp[2 * (r-1) + i2]);

			tridag(C, D, C1, rr, Ccr, r + 1);
			for (j3 = 0; j3 < r + 1; j3++)
			{
				Cc[x][j3] = Ccr[j3];
			}
		}

		//x=0
		for (j = 0; j < r + 1; j++)
		{
			Cavau[1] = (MTutemp[j] - Cmutemp[1][j]) / MTutemp[j];
			Cavad[0] = (MTdtemp[j] - Cmdtemp[0][j]) / MTdtemp[j];

			conu = -vm*Cavau[1] * Cmutemp[0][j] / dx[0];
			cond = vm*Cavad[0] * Cmdtemp[1][j] / dx[0];

			Cmu[0][j] = Cmutemp[0][j] + dt*(km*(MTutemp[j] - Cmutemp[0][j])*Cctemp[0][j] - kf*Cmutemp[0][j] + conu
				- Kau[j] * Cmutemp[0][j] - Cmutemp[0][j] / MTutemp[j] * k0*(MTutemp[j] - Cmutemp[x ][j]));
			Cmd[0][j] = Cmdtemp[0][j] + dt*(km*(MTdtemp[j] - Cmdtemp[0][j])*Cctemp[0][j] - kf*Cmdtemp[0][j] + cond
				- Kad[j] * Cmdtemp[0][j] - Cmdtemp[0][j] / MTdtemp[j] * k0*(MTdtemp[j] - Cmdtemp[0][j]));
		}
		for (i2 = 1; i2 < x; i2++)
		{
			for (j2 = 0; j2 < r + 1; j2++)
			{
				Cavau[i2] = (MTutemp[j2] - Cmutemp[i2][j2]) / MTutemp[j2];
				Cavau[i2 + 1] = (MTutemp[j2] - Cmutemp[i2 + 1][j2]) / MTutemp[j2];
				Cavad[i2] = (MTdtemp[j2] - Cmdtemp[i2][j2]) / MTdtemp[j2];
				Cavad[i2 - 1] = (MTdtemp[j2] - Cmdtemp[i2 - 1][j2]) / MTdtemp[j2];
				conu = vm*(Cavau[i2] * Cmutemp[i2 - 1][j2] - Cavau[i2 + 1] * Cmutemp[i2][j2]) / (dx[i2 - 1] + dx[i2]);
				cond = -vm*(Cavad[i2 - 1] * Cmdtemp[i2][j2] - Cavad[i2] * Cmdtemp[i2 + 1][j2]) / (dx[i2] + dx[i2 - 1]);
				Cmu[i2][j2] = Cmutemp[i2][j2] + dt*(km*(MTutemp[j2] - Cmutemp[i2][j2]) * Cctemp[i2][j2] - kf*Cmutemp[i2][j2] + conu
					- Kau[j2] * Cmutemp[i2][j2] - Cmutemp[i2][j2] / MTutemp[j2] * k0*(MTutemp[j2] - Cmutemp[x][j2]));
				Cmd[i2][j2] = Cmdtemp[i2][j2] + dt* (km*(MTdtemp[j2] - Cmdtemp[i2][j2])* Cctemp[i2][j2] - kf*Cmdtemp[i2][j2] + cond
					- Kad[j2] * Cmdtemp[i2][j2] - Cmdtemp[i2][j2] / MTdtemp[j2] * k0*(MTdtemp[j2] - Cmdtemp[0][j2]));
			}
		}
		for (j3 = 0; j3 < r + 1; j3++)//x
		{
			Cavau[x] = (MTutemp[j3] - Cmutemp[x][j3]) / MTutemp[j3];
			Cavad[x - 1] = (MTdtemp[j3] - Cmdtemp[x - 1][j3]) / MTdtemp[j3];
			conu = vm*Cavau[x] * Cmutemp[x - 1][j3] / dx[x - 1];
			cond = -vm*Cavad[x - 1] * Cmdtemp[x][j3] / dx[x - 1];
			Cmu[x][j3] = Cmutemp[x][j3] + dt*(km*(MTutemp[j3] - Cmutemp[x][j3]) * Cctemp[x][j3] - kf*Cmutemp[x][j3] + conu
				- Kau[j3] * Cmutemp[x][j3] - Cmutemp[x][j3] / MTutemp[j3] * k0*(MTutemp[j3] - Cmutemp[x][j3]));

			Cmd[x][j3] = Cmdtemp[x][j3] + dt*(km*(MTdtemp[j3] - Cmdtemp[x][j3]) * Cctemp[x][j3] - kf*Cmdtemp[x][j3] + cond
				- Kad[j3] * Cmdtemp[x][j3] - Cmdtemp[x][j3] / MTdtemp[j3] * k0*(MTdtemp[j3] - Cmdtemp[0][j3]));

		}



		for (j = 0; j < r + 1; j++)
		{

			Kmtu[j] = FKm(MTutemp[j],kp1,KMT);
			Kmtd[j] = FKm(MTdtemp[j], kp1, KMT);
			dMTu[j] = MT0temp[j] * (kp + Kmtu[j])*MT0temp[j] - k0*(MTutemp[j] - Cmutemp[x][j]) - Kau[j] * MTutemp[j];
			MTu[j] = MTutemp[j] + dt*dMTu[j];//
			dMTd[j] = MT0temp[j] * (kp + Kmtd[j])*MT0temp[j] - k0*(MTdtemp[j] - Cmdtemp[0][j]) - Kad[j] * MTdtemp[j];
			MTd[j] = MTdtemp[j] + dt*dMTd[j];//
			dMT0[j] = -dMTu[j] - dMTd[j];
			MT0[j] = MT0temp[j] + dt*dMT0[j];//
		}
		if (k % 5000 == 0 && k*dt >= 7200)
		{
			MTsum = 0;
			MTus = 0;
			MTds = 0;
			MT0s = 0;
			Cau = 0;
			Cad = 0;

			for (j = 0; j < r; j++)
			{
				
				MTus = MTus + 0.5*(MTu[j] + MTu[j + 1]) * 2 * pi*dr[j] * (rn[j] + 0.5*dr[j]);
				MTds = MTds + 0.5*(MTd[j] + MTd[j + 1]) * 2 * pi*dr[j] * (rn[j] + 0.5*dr[j]);
				MT0s = MT0s + 0.5*(MT0[j] + MT0[j + 1]) * 2 * pi*dr[j] * (rn[j] + 0.5*dr[j]);
				
			}
			MTsum = MTus + MTds + MT0s;

			fprintf(outMTu, "%lf\t", MTus);
			fprintf(outMTd, "%lf\t", MTds);
			fprintf(outMT0, "%lf\t", MT0s);
			fprintf(outMT, "%lf\t", MTsum);

			Csum = 0;
			Cmdsum = 0;
			Cmusum = 0;
			Ccsum1 = 0;


			for (i = 0; i < x + 1; i++)
			{
				Ccsum[i] = 0;
			}
				//kkk = kkk + 1;
				//sprintf(buffer, "Z%d.dat", kkk);
				//ofstream xxfile(buffer, ios::out);
				

			for (i = 0; i < x; i++)
			{
				for (j = 0; j < r; j++)
				{
					Ccsum[i] = Ccsum[i] + 0.25*(Cc[i][j] + Cc[i + 1][j] + Cc[i][j + 1] + Cc[i + 1][j + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] * dx[i];
					Cmds = 0.25*(Cmd[i][j] + Cmd[i][j + 1] + Cmd[i + 1][j] + Cmd[i + 1][j + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] * dx[i];
					Cmdsum = Cmdsum + Cmds;
					Cmus = 0.25*(Cmu[i][j] + Cmu[i][j + 1] + Cmu[i + 1][j] + Cmu[i + 1][j + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] * dx[i];
					Cmusum = Cmusum + Cmus;
					Ccsum1 = Ccsum1 + 0.25*(Cc[i][j] + Cc[i + 1][j] + Cc[i][j + 1] + Cc[i + 1][j + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j] * dx[i];
				}
				//xxfile << Ccsum[i] << "\t";
				Csum = Csum + Ccsum[i];
				fprintf(kemoCmd, "%lf\t", Cmds);
				fprintf(kemoCmu, "%lf\t", Cmus);
				fprintf(kemoCc, "%lf\t", Ccs);
			}//x=x+1




			fprintf(kemoCmd, "\n");
			fprintf(kemoCmu, "\n");
			fprintf(kemoCc, "\n");

			fprintf(outCmu, "%lf\t", Cmusum);
			fprintf(outCmd, "%lf\t", Cmdsum);

			Casum = 0;

			for (j = 0; j < r - 1; j++)
			{
				Casum = Casum + 0.5*(Ca[2 * j] + Ca[2 * (j + 1)] + Ca[2 * j + 1] + Ca[2 * (j + 1) + 1])
					* 2 * pi*(rn[j] + dr[j] / 2)*dr[j];
				Cad = Cad + 0.5*(Ca[2 * j] + Ca[2 * (j + 1)]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j];
				Cau = Cau + 0.5*(Ca[2 * j + 1] + Ca[2 * (j + 1) + 1]) * 2 * pi*(rn[j] + dr[j] / 2)*dr[j];

			}
			Casum = Casum + 0.5*(Ca[2 * (r - 1)] + Ca[2 * r]) * 2 * pi*(rn[j] + dr[r - 1] / 2)*dr[r - 1];
			Casum = Casum + 0.5*(Ca[2 * r - 1] + Ca[2 * r + x]) * 2 * pi*(rn[j] + dr[r - 1] / 2)*dr[r - 1];

			Cad = Cad + 0.5*(Ca[2 * (r - 1)] + Ca[2 * r]) * 2 * pi*(rn[j] + dr[r - 1] / 2)*dr[r - 1];
			Cau = Cau + 0.5*(Ca[2 * r - 1] + Ca[2 * r + x]) * 2 * pi*(rn[j] + dr[r - 1] / 2)*dr[r - 1];
			//x line
			for (i = 0; i < x; i++)
			{
				fprintf(kemoCa, "%lf\t", Ca[2 * r + i]);
				Casum = Casum + 0.5*(Ca[2 * r + i] + Ca[2 * r + i + 1]) * 2 * pi*rad*dx[i];
			}
			fprintf(kemoCa, "%lf\t", Ca[2 * r + x]);
			fprintf(kemoCa, "\n");
			Csum = Csum + Cmusum + Cmdsum + Casum;

			fprintf(outCau, "%lf\t", Cau);
			fprintf(outCad, "%lf\t", Cad);
			fprintf(outCa, "%lf\t", Casum);
			fprintf(out, "%lf\t", Csum);
			fprintf(outCc, "%lf\t", Ccsum1);

			fprintf(out, "\n");
	//		cout << k*dt << " Cc " << Ccsum1 << " Ca " << Casum << " Cmu " << Cmusum << " Cmd " << Cmdsum << " sum " << Casum + Ccsum1 + Cmusum + Cmdsum <<" Cau "<<Ca[2*4+1]<<" Cad "<<Ca[2*4]<< endl;
		//	cout<< " MTu " << MTu[4] << " MTd " << MTd[4] << " MT0 "<<MT0[4]<<" difu "<<MTu[4]-Cmu[4][4]<<" difd "<<MTd[4]-Cmd[4][4]<<endl;
		//	cout << k*dt << "  MTu  " << MTus << "  MTd  " << MTds << "  MT0  " << MT0s << "  sum  " << MTsum << endl;
			fprintf(outCar, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Ca[0], Ca[1], Ca[10], Ca[11], Ca[30], Ca[31],  MTd[0], MTu[0], MTd[10], MTu[10], MTd[30], MTu[30]);
	//		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\%lf\t%lf\t\%lf", Dc, vm, km, kf, ka, ka1, kd, kp, kp1, KMT, b, k0, kc, alpha, Kca);
			fprintf(time, "%lf \t", k*dt);
			fprintf(testratesca, "%lf\t%lf\n", FKa(Ca[0], kc, Kca), FKa(Ca[1], kc, Kca));
			fprintf(testratesMT, "%lf\t%lf\n", FKm(MTd[0], kp1, KMT), FKm(MTu[0], kp1,KMT));

		}
	}
	//FILE * outCaend;
	//outCaend = fopen("Ca1.dat", "w");
	//for (q = 0; q < ai; q++)
	//{
	//	fprintf(outCaend, "%lf\t", Ca[q]);
	//}
	fprintf(para, "Dc %.0f vm%.1f km%.1f kf%.1f ka%.1f ka1%.0f kd%.1f kp%.1f kp1%.0f KMT%.0f b%.0f k0%.3f kc%.0f a%.0f kca%.2f L%.0f r%.2f", Dc, vm, km, kf, ka, ka1, kd, kp, kp1, KMT, b, k0, kc, alpha, Kca,L,rad);
	fclose(para);
	fprintf(pa, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\%lf\t%lf\t\%lf\f%lf\t%lf", Dc, vm, km, kf, ka, ka1, kd, kp, kp1, KMT, b, k0, kc, alpha, Kca,L,rad);
	fclose(pa);
	fclose(time);
	fclose(kemoCmu);
	fclose(kemoCmd);
	fclose(outCau);
	fclose(outCad);
	fclose(outCc);
	fclose(outMTu);
	fclose(outMTd);
	fclose(outMT0);
	fclose(outMT);
	fclose(outCmu);
	fclose(outCmd);
	fclose(outCa);
	fclose(kemoCa);
	fclose(outCar);

	fclose(testratesca);
	fclose(testratesMT);
	fclose(kemoCc);
	final = clock() - init;
	cout << (double)final / ((double)CLOCKS_PER_SEC);

	//cout << "set1";
	return 0;
}

