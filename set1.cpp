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

