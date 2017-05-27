
#include <iostream>					// for handeling IO stuff
#include <cmath>					// for general purpose math
#include <armadillo>				// for minear algebra
#include <complex>					// for complex numvers
#include <chrono>					// for measuring time

double pi = 4*atan(1);				// pi
double umass = 1.660539040e-27;		// atomic mass unit
double spring = 572.936;		// spring constant for molecular hydrogen
double hmass = 1.00794*umass;		// mass of hydrogen
double cutoff = 50e-12;				// cut off distance for potential
double V0 = 0;						// potentil at bottom of well
double hbar = 1.0545718e-34;
double omega = sqrt(spring/hmass);
double sqrtmohbar = sqrt(hmass*omega/hbar);

using namespace std;				// lazyness
using namespace arma;				// lazyness again

double V(double pos)				// calculate harmonic potential
{
	return 0.5*spring*pow(pos,2);	
}

double Vshift(double pos, double C1, double C2)
{
	return C1*exp(-C2*pow(pos*sqrtmohbar,2));
}

double aneigenval(int n)
{
	return hbar*omega*(n+0.5);
}

void hamiltonian(mat* H, double stepsize, int steps)		// constructing hamiltonian
{
	double c0 = -pow(hbar,2)/(2*hmass)*30.0/(-12.0*pow(stepsize,2));		// defining the c's
	double c1 = -pow(hbar,2)/(2*hmass)*(-16.0)/(-12.0*pow(stepsize,2));
	double c2 = -pow(hbar,2)/(2*hmass)*1.0/(-12.0*pow(stepsize,2));

	int m = 1; int n = 2;
	for(int k = 0; k < steps; k++)					// constructing matrix elements
	{
		//cout << -cutoff+k*stepsize << endl;
		(*H)(k,k) = c0 + V(-cutoff+k*stepsize);
		(*H)(m-1,m) = c1;
		(*H)(m,m-1) = c1;
		(*H)(n-2,n) = c2;
		(*H)(n,n-2) = c2;
		
		if (m < steps-1) m++;
		if (n < steps-1) n++;
	}
	//(*H)(steps-2,steps-1) = c1;
	//(*H)(steps-1,steps-2) = c1;
	//(*H)(steps-1,steps-1) = c0;

	//mat L, U, P;
	
	//lu(L,U,P,(*H));

	//cout << "Hamiltonian: " << endl;
	//(*H).print();

	//cout << "Upper triagonal matrix: " << endl;
	//U.print();

	vec eigenvals;
	mat eigenvecs;

	eig_sym(eigenvals,eigenvecs,(*H));

	//cout << "Eigen values: " << endl;
	//eigenvals.print();
	//cout << "Eigen vectors: " << endl;
	//eigenvecs.print();

	eigenvecs.save("eigenvectors",raw_ascii);
	eigenvals.save("eigenvalues",raw_ascii);

	//vec a = eigenvecs.col(0);
	//vec b = eigenvecs.col(1);
	
	//cout << "orthogonality: " << dot(a,b) << endl;

	//mat sanity1 = (*H)*eigenvecs;
	//mat sanity2 = eigenvals*eigenvecs;

	//mat diff = sanity1 - sanity2;

	//mat eigendiag = diagmat(eigenvals);
	//mat inveigenvecs = inv(eigenvecs);

	//eigendiag.print();

	//mat A = eigenvecs*eigendiag*inveigenvecs;
	
	//mat L2, U2, P2;
	//lu(L2,U2,P2,A);

	//cout << "Should be same as triangular above: " << endl;
	//U2.print();

	//mat diff = A-(*H);

	//cout << "Should be Hamiltonian too:: " << endl;
	//A.print();
	//diff.print();
	
	
}

void Hshifted(mat* H, double stepsize, int steps, double C1, double C2)		// constructing hamiltonian
{
	double c0 = -pow(hbar,2)/(2*hmass)*30.0/(-12.0*pow(stepsize,2));		// defining the c's
	double c1 = -pow(hbar,2)/(2*hmass)*(-16.0)/(-12.0*pow(stepsize,2));
	double c2 = -pow(hbar,2)/(2*hmass)*1.0/(-12.0*pow(stepsize,2));

	int m = 1; int n = 2;
	for(int k = 0; k < steps; k++)					// constructing matrix elements
	{
		(*H)(k,k) = c0 + V(-cutoff+k*stepsize) + Vshift(-cutoff+k*stepsize,C1,C2);
		(*H)(m-1,m) = c1;
		(*H)(m,m-1) = c1;
		(*H)(n-2,n) = c2;
		(*H)(n,n-2) = c2;
		
		if (m < steps-1) m++;
		if (n < steps-1) n++;
	}


	//cout << "Hamiltonian: " << endl;
	//(*H).print();

	//cout << "Upper triagonal matrix: " << endl;
	//U.print();

	vec eigenvals;
	mat eigenvecs;

	eig_sym(eigenvals,eigenvecs,(*H));

	//cout << "Eigen values: " << endl;
	//eigenvals.print();
	//cout << "Eigen vectors: " << endl;
	//eigenvecs.print();

	eigenvecs.save("eigenvectors2",raw_ascii);
}

void invpowit(mat* H, vec* esteigen, double pres, double diff)
{
	//cout << "input " << endl;
	//(*esteigen).print();
	
	double acc1 = accu((*esteigen));

	(*esteigen) = solve((*H),(*esteigen));
	(*esteigen) = normalise((*esteigen));

	//cout << "output " << endl;
	//(*esteigen).print();
	
	double acc2 = accu((*esteigen));
	double diff2 = abs(acc1 + acc2);

	//cout << diff << " " << diff2 << " " << diff-diff2 << endl;

	if (abs(diff - diff2) < 1e14)
	{diff = diff2 - 1e25;}
	else
	{diff = diff2;}

	//cout << "diff: " << diff << endl;
	//if (pres != 0)
	if (diff > pres)
		invpowit(H,esteigen,pres,diff);
}

void invpowitfast(mat* H, vec* esteigen, double pres)
{
	//cout << "input " << endl;
	//(*esteigen).print();

	double diff = 1; 
	double acc1; double acc2;
	while (diff > pres)
	{
		acc1 = accu((*esteigen));

		(*esteigen) = solve((*H),(*esteigen));
		(*esteigen) = normalise((*esteigen));

		//cout << "output " << endl;
		//(*esteigen).print();
	
		acc2 = accu((*esteigen));
		diff = abs(acc1 + acc2);

		//cout << "diff: " << diff << endl;
		//if (pres != 0)
	}
}

double shiftinvpowit(mat* H, double pres, double shift, int steps)
{
	vec esteigenvec(steps,fill::ones);
	vec shiftvec(steps,fill::zeros);
	shiftvec -= -shift;
	mat shiftmat = diagmat(shiftvec);
	(*H) -= shiftmat;
	invpowit(H,&esteigenvec,pres,0);


	mat invH = inv((*H));
	double esteigenval = as_scalar(esteigenvec.t()*invH*esteigenvec)/as_scalar(esteigenvec.t()*esteigenvec);
	
	esteigenval = 1.0/esteigenval + shift;

	return esteigenval;
}

int main()
{
	int steps = 500;
	double pres = 1e-15;
	double stepsize = 2*cutoff/steps;
	double shift = 6.5e-19;
	double shift1 = 1e-18;
	double shift2 = 1e0;

	mat H(steps,steps,fill::zeros);					// initialize matrix of zeroes
	
	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::    now(); 
	
	hamiltonian(&H,stepsize,steps);
	
	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::    now();

	mat eigenvals;
	eigenvals.load("eigenvalues");	

	double time1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 )    .count();


	mat Hshift(steps,steps,fill::zeros);					// initialize matrix of zeroes
	Hshifted(&Hshift,stepsize,steps,shift1,shift2);

	mat HU, HL, HP;
	lu(HU,HL,HP,H);
	mat H2 = HL*HU;
	
	chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::    now(); 

	double esteigenval = shiftinvpowit(&H2, pres, shift, steps);
	
	chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::    now(); 
	
	double time2 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 )    .count();
	
	int ni = 10;
	cout << "calculated eigenval: " << eigenvals(ni) << endl << "diff from analytic: " << abs(eigenvals(ni) - aneigenval(ni)) << endl;
	cout << "estimated eigenval: " << esteigenval << endl << "diff from analytic: " << abs(esteigenval-aneigenval(ni)) << endl;
	cout << "methods diff by: " << abs(esteigenval - eigenvals(ni)) << endl;
	
	cout << "time1 " << time1 << endl << "time2 " << time2 << endl;
	
	return 0;
}














