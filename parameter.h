const double pi = 3.14159265358979323846264338327950288;
//const int x = 140, r =30;
const int L = 1;//Keep
const int L1 =127;
const double rad = sqrt(10 * 10 * 130 / L1);//keep
const int r1 = 10;
const int r2 = 10;
const int r3 = 10;
const int r4 = 10;
const int r5 = 10;

const double dr1 = rad / 50;
const double dr2 = rad / 50;
const double dr3 = rad / 50;
const double dr4 = rad / 50;
const double dr5 = rad / 50;



const double dx1 = 0.01;
const double dx2 = 0.02;
const double dx3 = 0.05;
const double dx4 = 0.2;
const double dx5 = 1;
const int x1 = 10;//first and last 10
const int x2 = 5;// 5 later
const int x3 = 4;//4 later
const int x4 = 3;//3 later
const int x5 = L1 - 2;//
const int x = 2 * (x1 + x2 + x3 + x4) + x5, r = r1 + r2 + r3 + r4 + r5;
//const double dr1 = 0.05, r1 = 4;//last 10
//const double dr2 = 0.05, r2 = 4;// 5 later
//const double dr3 = 0.2, r3 = 3;//4 later
//const double dr4 = 0.5, r4 = 2;//3 later
//const double dr5 = 1, r5 =3;//4 other
//double xcheck = 2 * (x1 + x2 + x3 + x4) + x5,rcheck= r1+r2+r3+r4+r5;
const int ai = 2 * r + x + 1;
const double Ccall = 1200;

const double tf = 36000;


const double dt = 0.001;//I would keep it that way for now...
const double Dc0 = 10;// 1-10 typical Cell Signalling : Changing Shape Changes the Signal
const double vm0 = 0.3;// Tracking Single Kinesin Molecules in the Cytoplasm of Mammalian Cells
const double km0 = 1;
const double kf0 = 1 / 1.5;//Tracking Single Kinesin Molecules in the Cytoplasm of Mammalian Cells
const double ka0 = 0.3;//
const double ka10 = 8;//
const double kd0 = 0.3;// 0.1-100.... Cell Signalling : Changing Shape Changes the Signal
const double kp0 = 0.1;//review paper
const double kp10 = 3;
const double KMT0 = 6;
const double b = 2;
const double k00 = 0.017;//old paper
const double kc0 = 400;//review
const double alpha = 10;
const double Kca0 = 0.3;


//const double vm = 0;
//const double km = 0;
//const double kf = 0;
//const double ka = 0;
//const double ka1 = 0;
//const double kd = 0;
//const double kp = 0;
//const double kp1 = 0;
//const double KMT = 0;
//const double b = 0;
//const double k0 = 0;
//const double kc = 0;
//const double alpha = 0;// 5;
//const double Kca = 0;// 0.5;