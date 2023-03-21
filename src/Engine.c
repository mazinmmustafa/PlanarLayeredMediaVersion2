//
#include "Engine.h"

complex double sqrtm(complex double z){
	complex double w=csqrt(z);
	if (cimag(w)<=0.0){
		return +w;
	}else{
		return -w;
	}
}

void Parameters_k_rho(int n, char p, complex double k_rho, complex double *kz, 
	complex double *Z, complex double *Theta, MLayers *myConfig){
	assert(myConfig->Layers!=NULL);
	int N=myConfig->N;
	assert(n<=N&&n>0);
	assert(p=='e'||p=='h');
	double omega=2.0*pi*myConfig->freq;
	double d=myConfig->Layers[n-1].d;
	complex double eps=myConfig->Layers[n-1].eps;
	complex double mu=myConfig->Layers[n-1].mu;
	complex double k=myConfig->Layers[n-1].k;
	*(kz) = sqrtm(k*k-k_rho*k_rho);
	*(Theta) = (*kz)*d;
	if (p=='e'){
		*(Z) = (*kz)/(omega*eps0*eps);
	}
	if (p=='h'){
		*(Z) = (omega*mu0*mu)/(*kz);
	}
}

// complex double Gamma_d(int n, char p, complex double k_rho, MLayers *myConfig){
// 	assert(myConfig->Layers!=NULL);
// 	complex double j=csqrt(-1.0);
// 	int N=myConfig->N;
// 	assert(n<=N&&n>0);
// 	if (n==N){
// 		return myConfig->Gamma_d;
// 	}else{
// 		complex double kz_n, kz_m;
// 		complex double Z_n, Z_m;
// 		complex double Theta_n, Theta_m;
// 		complex double sigma_s=myConfig->Layers[n-1].sigma_s;
// 		complex double Gamma, Omega;
// 		Parameters_k_rho(n+0, p, k_rho, &kz_n, &Z_n, &Theta_n, myConfig);
// 		Parameters_k_rho(n+1, p, k_rho, &kz_m, &Z_m, &Theta_m, myConfig);
// 		Gamma = (Z_m-Z_n)/(Z_m+Z_n);
// 		Omega = (Z_m*Z_n)/(Z_m+Z_n);
// 		complex double A, B, C, D;
// 		A = Gamma-Omega*sigma_s;
// 		B = 1.0-Omega*sigma_s;
// 		C = 1.0+Omega*sigma_s;
// 		D = Gamma+Omega*sigma_s;
// 		complex double Exp=cexp(-j*2.0*Theta_m)*Gamma_d(n+1, p, k_rho, myConfig);
// 		return (A+B*Exp)/(C+D*Exp);
// 	}
// }

// complex double Gamma_u(int n, char p, complex double k_rho, MLayers *myConfig){
// 	assert(myConfig->Layers!=NULL);
// 	complex double j=csqrt(-1.0);
// 	int N=myConfig->N;
// 	assert(n<=N&&n>0);
// 	if (n==1){
// 		return myConfig->Gamma_u;
// 	}else{
// 		complex double kz_n, kz_m;
// 		complex double Z_n, Z_m;
// 		complex double Theta_n, Theta_m;
// 		complex double sigma_s=myConfig->Layers[n-2].sigma_s;
// 		complex double Gamma, Omega;
// 		Parameters_k_rho(n+0, p, k_rho, &kz_n, &Z_n, &Theta_n, myConfig);
// 		Parameters_k_rho(n-1, p, k_rho, &kz_m, &Z_m, &Theta_m, myConfig);
// 		Gamma = (Z_m-Z_n)/(Z_m+Z_n);
// 		Omega = (Z_m*Z_n)/(Z_m+Z_n);
// 		complex double A, B, C, D;
// 		A = Gamma-Omega*sigma_s;
// 		B = 1.0-Omega*sigma_s;
// 		C = 1.0+Omega*sigma_s;
// 		D = Gamma+Omega*sigma_s;
// 		complex double Exp=cexp(-j*2.0*Theta_m)*Gamma_u(n-1, p, k_rho, myConfig);
// 		return (A+B*Exp)/(C+D*Exp);
// 	}
// }

complex double Gamma_d(int n, char p, complex double k_rho, MLayers *myConfig){
	assert(myConfig->Layers!=NULL);
	complex double j=csqrt(-1.0);
	int N=myConfig->N;
	assert(n<=N&&n>0);
	complex double kz_n, kz_m;
	complex double Z_n, Z_m;
	complex double Theta_n, Theta_m;
	complex double sigma_s;
	complex double Omega;
	complex double A, B, C, D;
	complex double Gamma, Exp;
	complex double Refl=myConfig->Gamma_d;
	for (int i=(N-1); i>=n; i--){
		sigma_s=myConfig->Layers[i-2].sigma_s;
		Parameters_k_rho(i+0, p, k_rho, &kz_n, &Z_n, &Theta_n, myConfig);
		Parameters_k_rho(i+1, p, k_rho, &kz_m, &Z_m, &Theta_m, myConfig);
		Gamma = (Z_m-Z_n)/(Z_m+Z_n);
		Omega = (Z_m*Z_n)/(Z_m+Z_n);
		A = Gamma-Omega*sigma_s;
		B = 1.0-Omega*sigma_s;
		C = 1.0+Omega*sigma_s;
		D = Gamma+Omega*sigma_s;
		Exp = Refl*cexp(-j*2.0*Theta_m);
		Refl = (A+B*Exp)/(C+D*Exp);
	}
	return Refl;
}

complex double Gamma_u(int n, char p, complex double k_rho, MLayers *myConfig){
	assert(myConfig->Layers!=NULL);
	complex double j=csqrt(-1.0);
	int N=myConfig->N;
	assert(n<=N&&n>0);
	complex double kz_n, kz_m;
	complex double Z_n, Z_m;
	complex double Theta_n, Theta_m;
	complex double sigma_s;
	complex double Omega;
	complex double A, B, C, D;
	complex double Gamma, Exp;
	complex double Refl=myConfig->Gamma_u;
	for (int i=2; i<=n; i++){
		sigma_s=myConfig->Layers[i-2].sigma_s;
		Parameters_k_rho(i+0, p, k_rho, &kz_n, &Z_n, &Theta_n, myConfig);
		Parameters_k_rho(i-1, p, k_rho, &kz_m, &Z_m, &Theta_m, myConfig);
		Gamma = (Z_m-Z_n)/(Z_m+Z_n);
		Omega = (Z_m*Z_n)/(Z_m+Z_n);
		A = Gamma-Omega*sigma_s;
		B = 1.0-Omega*sigma_s;
		C = 1.0+Omega*sigma_s;
		D = Gamma+Omega*sigma_s;
		Exp = Refl*cexp(-j*2.0*Theta_m);
		Refl = (A+B*Exp)/(C+D*Exp);
	}
	return Refl;
}

int selectLayer(double z, MLayers *myConfig){
	int N=myConfig->N;
	int n=N;
	for (int i=N-1; i>0; i--){
		if (z>=myConfig->Layers[i-1].zd){
			n = i;
		}
	}
	return n;
}

void TLGF(double z, double z_, complex double k_rho, char p, char s, int option, 
	MLayers *myConfig, complex double *V, complex double *I){
	assert(myConfig->Layers!=NULL);
	assert((s=='v')||(s='i'));
	assert((option==0)||(option==1));
	complex double j=csqrt(-1.0);
	(*V) = 0.0;
	(*I) = 0.0;
	int n=selectLayer(z, myConfig);
	int m=selectLayer(z_, myConfig);
	complex double v, i;
	if (s=='v'){
		v = 1.0;
		i = 0.0;
	}
	if (s=='i'){
		v = 0.0;
		i = 1.0;
	}
	complex double Vm_p, Vm_m;
	double zu, zd;
	zu = myConfig->Layers[m-1].zu;
	zd = myConfig->Layers[m-1].zd;
	complex double Theta_m;
	complex double Gamma_m_u, Gamma_m_d;
	Gamma_m_d = Gamma_d(m, p, k_rho, myConfig);
	Gamma_m_u = Gamma_u(m, p, k_rho, myConfig);
	complex double Dm, Zm, kz_m;
	Parameters_k_rho(m, p, k_rho, &kz_m, &Zm, &Theta_m, myConfig);
	Dm = 1.0-Gamma_m_d*Gamma_m_u*cexp(-j*2*Theta_m);
	Vm_p = (1.0-Gamma_m_d*cexp(-j*2.0*kz_m*(z_-zd)))*v;
	Vm_p+= (1.0+Gamma_m_d*cexp(-j*2.0*kz_m*(z_-zd)))*i*Zm;
	Vm_p/= 2.0*Dm;
	Vm_m = (1.0+Gamma_m_u*cexp(-j*2.0*kz_m*(zu-z_)))*i*Zm;
	Vm_m-= (1.0-Gamma_m_u*cexp(-j*2.0*kz_m*(zu-z_)))*v;
	Vm_m/= 2.0*Dm;
	complex double Vm0_p, Vm0_m;
	Vm0_p = 0.5*(v+i*Zm);
	Vm0_m = 0.5*(i*Zm-v);
	if (m==n){
		complex double Gamma_m_u_, Gamma_m_d_;
		Gamma_m_u_ = Gamma_m_u*cexp(-j*2.0*kz_m*(zu-z_));
		Gamma_m_d_ = Gamma_m_d*cexp(-j*2.0*kz_m*(z_-zd));
		if (z>=z_){
			if (option==0){
				(*V) = +Vm_p*(cexp(-j*kz_m*(z-z_))+Gamma_m_u_*cexp(+j*kz_m*(z-z_)));
				(*I) = +Vm_p*(cexp(-j*kz_m*(z-z_))-Gamma_m_u_*cexp(+j*kz_m*(z-z_)))/Zm;
			}
			if (option==1){
				(*V) = (Vm_p-Vm0_p)*cexp(-j*kz_m*(z-z_))+Vm_p*Gamma_m_u_*cexp(+j*kz_m*(z-z_));
				(*I) = ((Vm_p-Vm0_p)*cexp(-j*kz_m*(z-z_))-Vm_p*Gamma_m_u_*cexp(+j*kz_m*(z-z_)))/Zm;
			}
		}
		if (z<z_){
			if (option==0){
				(*V) = +Vm_m*(cexp(+j*kz_m*(z-z_))+Gamma_m_d_*cexp(-j*kz_m*(z-z_)));
				(*I) = -Vm_m*(cexp(+j*kz_m*(z-z_))-Gamma_m_d_*cexp(-j*kz_m*(z-z_)))/Zm;
			}
			if (option==1){
				(*V) = (Vm_m-Vm0_m)*cexp(+j*kz_m*(z-z_))+Vm_m*Gamma_m_d_*cexp(-j*kz_m*(z-z_));
				(*I) = ((Vm0_m-Vm_m)*cexp(+j*kz_m*(z-z_))+Vm_m*Gamma_m_d_*cexp(-j*kz_m*(z-z_)))/Zm;
			}
		}
	}
	if (m<n){
		complex double Gamma_n_d;
		complex double Vn_m, kz_n, Zn, Theta_n;
		Gamma_n_d = Gamma_d(m+1, p, k_rho, myConfig);
		Parameters_k_rho(m+1, p, k_rho, &kz_n, &Zn, &Theta_n, myConfig);
		Vn_m = Vm_m;
		Vn_m*=(1.0+Gamma_m_d)*cexp(-j*kz_m*(z_-zd));
		Vn_m/=(1.0+Gamma_n_d*cexp(-j*2.0*Theta_n));
		for (int i=m+2; i<=n; i++){
			Gamma_m_d = Gamma_n_d;
			Gamma_n_d = Gamma_d(i, p, k_rho, myConfig);
			kz_m = kz_n; 
			Theta_m = Theta_n;
			Parameters_k_rho(i, p, k_rho, &kz_n, &Zn, &Theta_n, myConfig);
			Vn_m*=(1.0+Gamma_m_d)*cexp(-j*Theta_m);
			Vn_m/=(1.0+Gamma_n_d*cexp(-j*2.0*Theta_n));
		}
		zd = myConfig->Layers[n-1].zd;
		(*V) = +Vn_m*cexp(-j*Theta_n)*(cexp(+j*kz_n*(z-zd))+Gamma_n_d*cexp(-j*kz_n*(z-zd)));
		(*I) = -Vn_m*cexp(-j*Theta_n)*(cexp(+j*kz_n*(z-zd))-Gamma_n_d*cexp(-j*kz_n*(z-zd)))/Zn;
	}
	if (m>n){
		complex double Gamma_n_u;
		complex double Vn_p, kz_n, Zn, Theta_n;
		Gamma_n_u = Gamma_u(m-1, p, k_rho, myConfig);
		Parameters_k_rho(m-1, p, k_rho, &kz_n, &Zn, &Theta_n, myConfig);
		Vn_p = Vm_p;
		Vn_p*=(1.0+Gamma_m_u)*cexp(-j*kz_m*(zu-z_));
		Vn_p/=(1.0+Gamma_n_u*cexp(-j*2.0*Theta_n));
		for (int i=m-2; i>=n; i--){
			Gamma_m_u = Gamma_n_u;
			Gamma_n_u = Gamma_u(i, p, k_rho, myConfig);
			kz_m = kz_n; 
			Theta_m = Theta_n;
			Parameters_k_rho(i, p, k_rho, &kz_n, &Zn, &Theta_n, myConfig);
			Vn_p*=(1.0+Gamma_m_u)*cexp(-j*Theta_m);
			Vn_p/=(1.0+Gamma_n_u*cexp(-j*2.0*Theta_n));
		}
		zu = myConfig->Layers[n-1].zu;
		(*V) = +Vn_p*cexp(-j*Theta_n)*(cexp(-j*kz_n*(z-zu))+Gamma_n_u*cexp(+j*kz_n*(z-zu)));
		(*I) = +Vn_p*cexp(-j*Theta_n)*(cexp(-j*kz_n*(z-zu))-Gamma_n_u*cexp(+j*kz_n*(z-zu)))/Zn;
	}
}

double findMax(MLayers *myConfig){
	double max=0.0;
	double temp=0.0;
	for (int n=0; n<myConfig->N; n++){
		temp = cabs(myConfig->Layers[n].k);
		if (temp>max){
			max = temp;
		}
	}
	return (max+myConfig->k0);
}

void detour_k_rho(complex double t, double a, double b, 
	complex double *k_rho, complex double *d_k_rho){
	assert(a!=0.0);
	complex double j=csqrt(-1.0);
	(*k_rho) = t+j*b*csin(pi*t/a);
	(*d_k_rho) = 1.0+j*(pi*b/a)*ccos(pi*t/a);
}

double findDetourParameters(double rho, double z, double z_, MLayers *myConfig){
	double k0=myConfig->k0;
	double b;
	if (rho>fabs(z-z_)){
		b = k0 < (1.0/rho) ? k0 : (1.0/rho);
	}else{
		b = k0;
	}
	return b;
}

complex double SI(complex double G(complex double, void*), void *args, 
	double a, double tol, int kmax){
	int flag;
	return QuadLAdaptive1D_32(G, args, 0.0, a, tol, kmax, 0, 0.0, &flag);
}

PlaneWaveFields PlaneWave(complex double ETheta, complex double EPhi,
	double theta_i, double phi_i, double x, double y, double z, MLayers *myConfig){
	assert(myConfig->Layers!=NULL);
	double z0=myConfig->Layers[0].zu;
	assert(cimag(myConfig->Layers[0].k)==0.0);
	double k1=creal(myConfig->Layers[0].k);
	double k_rho=k1*sin(theta_i);
	double ki_rho_dot=k_rho*cos(phi_i)*x+k_rho*sin(phi_i)*y;
	complex double Ze1=myConfig->Layers[0].eta*cos(theta_i);
	complex double Zh1=myConfig->Layers[0].eta/cos(theta_i);
	complex double eps, mu, eps1, mu1;
	int n=selectLayer(z, myConfig);
	eps1 = myConfig->Layers[0].eps;
	eps = myConfig->Layers[n-1].eps;
	mu1 = myConfig->Layers[0].mu;
	mu = myConfig->Layers[n-1].mu;
	complex double j=csqrt(-1.0);
	complex double Exp=-2.0*cexp(+j*ki_rho_dot)*cexp(+j*k1*cos(theta_i)*z0);
	complex double V_e_v, V_h_v, I_e_v, I_h_v;
	TLGF(z, z0, k_rho, 'e', 'v', 0, myConfig, &V_e_v, &I_e_v);
	TLGF(z, z0, k_rho, 'h', 'v', 0, myConfig, &V_h_v, &I_h_v);	
	complex double Ex, Ey, Ez;
	Ex =+ETheta*V_e_v*cos(theta_i)*(cos(phi_i))
		+EPhi*V_h_v*(-sin(phi_i));
	Ey =+ETheta*V_e_v*cos(theta_i)*(sin(phi_i))
		+EPhi*V_h_v*(+cos(phi_i));
	Ez =-ETheta*I_e_v*Ze1*(eps1/eps)*sin(theta_i);
	Field E={Ex*Exp, Ey*Exp, Ez*Exp};
	complex double Hx, Hy, Hz;
	Hx =+ETheta*I_e_v*cos(theta_i)*(-sin(phi_i))
		-EPhi*I_h_v*(cos(phi_i));
	Hy =+ETheta*I_e_v*cos(theta_i)*(+cos(phi_i))
		-EPhi*I_h_v*(sin(phi_i));
	Hz =+EPhi*V_h_v*(1.0/Zh1)*(mu1/mu)*tan(theta_i);
	Field H={Hx*Exp, Hy*Exp, Hz*Exp};
	PlaneWaveFields results={E, H};
	return results;
}

FarField FarFieldDipoleJ(Dipole J, double theta_s, double phi_s, double r, MLayers *myConfig){
	assert(myConfig->Layers!=NULL);
	int N=myConfig->N;
	double z0=myConfig->Layers[0].zu;
	double zN=myConfig->Layers[N-1].zd;
	double k1=creal(myConfig->Layers[0].k);
	double kN=creal(myConfig->Layers[N-1].k);
	double k1_rho=k1*sin(theta_s);
	double kN_rho=kN*sin(theta_s);
	double x_=J.x_;
	double y_=J.y_;
	double z_=J.z_;
	complex double Ilx=J.Ilx;
	complex double Ily=J.Ily;
	complex double Ilz=J.Ilz;
	int m=selectLayer(z_, myConfig);
	complex double eps_;
	eps_ = myConfig->Layers[m-1].eps;
	double k1s_rho_dot=k1_rho*cos(phi_s)*x_+k1_rho*sin(phi_s)*y_;
	double kNs_rho_dot=kN_rho*cos(phi_s)*x_+kN_rho*sin(phi_s)*y_;
	complex double eta1=myConfig->Layers[0].eta;
	complex double etaN=myConfig->Layers[N-1].eta;
	complex double eps1, epsN;
	eps1 = myConfig->Layers[0].eps;
	epsN = myConfig->Layers[N-1].eps;
	complex double j=csqrt(-1.0);
	complex double Exp1=(-j*k1*cexp(-j*k1*r)/(2.0*pi*r))*
						cexp(+j*k1s_rho_dot)*cexp(+j*k1*cos(theta_s)*z0);
	complex double ExpN=(-j*kN*cexp(-j*kN*r)/(2.0*pi*r))*
						cexp(+j*kNs_rho_dot)*cexp(+j*kN*cos(theta_s)*zN);
	complex double V_e_v, V_e_i, V_h_i, I;
	
	int flag=0;
	complex double E_theta, E_phi;
	FarField E;
	if ((theta_s>=0.0)&&(theta_s<=pi/2.0)){
		if (cimag(myConfig->Layers[0].k)<0.0){
			E.theta = 0.0;
			E.phi = 0.0;
		}else{
			TLGF(z0, z_, k1_rho, 'e', 'v', 0, myConfig, &V_e_v, &I);
			TLGF(z0, z_, k1_rho, 'e', 'i', 0, myConfig, &V_e_i, &I);	
			TLGF(z0, z_, k1_rho, 'h', 'i', 0, myConfig, &V_h_i, &I);
			E_theta =+V_e_i*(cos(phi_s)*Ilx+sin(phi_s)*Ily)
					 -V_e_v*eta1*(eps1/eps_)*sin(theta_s)*(Ilz);
			E_phi =+V_h_i*cos(theta_s)*(-sin(phi_s)*Ilx+cos(phi_s)*Ily);
			E.theta = E_theta*Exp1;
			E.phi = E_phi*Exp1;	
		}
		flag++;
	}
	if ((theta_s>pi/2.0)&&(theta_s<=pi)){
		if (cimag(myConfig->Layers[N-1].k)<0.0){
			E.theta = 0.0;
			E.phi = 0.0;
		}else{
			TLGF(zN, z_, kN_rho, 'e', 'v', 0, myConfig, &V_e_v, &I);
			TLGF(zN, z_, kN_rho, 'e', 'i', 0, myConfig, &V_e_i, &I);	
			TLGF(zN, z_, kN_rho, 'h', 'i', 0, myConfig, &V_h_i, &I);	
			E_theta =+V_e_i*(cos(phi_s)*Ilx+sin(phi_s)*Ily)
					 -V_e_v*etaN*(epsN/eps_)*sin(theta_s)*(Ilz);
			E_phi =+V_h_i*cos(theta_s)*(-sin(phi_s)*Ilx+cos(phi_s)*Ily);
			E.theta = E_theta*ExpN;
			E.phi = E_phi*ExpN;
		}
		flag++;
	}
	assert(flag==1);
	return E;
}

FarField FarFieldDipoleM(Dipole M, double theta_s, double phi_s, double r, MLayers *myConfig){
	assert(myConfig->Layers!=NULL);
	int N=myConfig->N;
	double z0=myConfig->Layers[0].zu;
	double zN=myConfig->Layers[N-1].zd;
	double k1=creal(myConfig->Layers[0].k);
	double kN=creal(myConfig->Layers[N-1].k);
	double k1_rho=k1*sin(theta_s);
	double kN_rho=kN*sin(theta_s);
	double x_=M.x_;
	double y_=M.y_;
	double z_=M.z_;
	complex double Ilx=M.Ilx;
	complex double Ily=M.Ily;
	complex double Ilz=M.Ilz;
	int m=selectLayer(z_, myConfig);
	complex double mu_;
	mu_ = myConfig->Layers[m-1].mu;
	double k1s_rho_dot=k1_rho*cos(phi_s)*x_+k1_rho*sin(phi_s)*y_;
	double kNs_rho_dot=kN_rho*cos(phi_s)*x_+kN_rho*sin(phi_s)*y_;
	complex double eta1=myConfig->Layers[0].eta;
	complex double etaN=myConfig->Layers[N-1].eta;
	complex double mu1, muN;
	mu1 = myConfig->Layers[0].mu;
	muN = myConfig->Layers[N-1].mu;
	complex double j=csqrt(-1.0);
	complex double Exp1=(+j*k1*cexp(-j*k1*r)/(2.0*pi*r))*
						cexp(+j*k1s_rho_dot)*cexp(+j*k1*cos(theta_s)*z0);
	complex double ExpN=(+j*kN*cexp(-j*kN*r)/(2.0*pi*r))*
						cexp(+j*kNs_rho_dot)*cexp(+j*kN*cos(theta_s)*zN);
	complex double V_e_v, V_h_i, V_h_v, I;
	
	int flag=0;
	complex double E_theta, E_phi;
	FarField E;
	if ((theta_s>=0.0)&&(theta_s<=pi/2.0)){
		if (cimag(myConfig->Layers[0].k)<0.0){
			E.theta = 0.0;
			E.phi = 0.0;
		}else{
			TLGF(z0, z_, k1_rho, 'e', 'v', 0, myConfig, &V_e_v, &I);
			TLGF(z0, z_, k1_rho, 'h', 'i', 0, myConfig, &V_h_i, &I);	
			TLGF(z0, z_, k1_rho, 'h', 'v', 0, myConfig, &V_h_v, &I);	
			E_theta =-V_e_v*(-sin(phi_s)*Ilx+cos(phi_s)*Ily);
			E_phi = +V_h_v*cos(theta_s)*(cos(phi_s)*Ilx+sin(phi_s)*Ily)
					-V_h_i*(1.0/eta1)*(mu1/mu_)*sin(theta_s)*cos(theta_s)*(Ilz);
			E.theta = E_theta*Exp1;
			E.phi = E_phi*Exp1;	
		}
		flag++;
	}
	if ((theta_s>pi/2.0)&&(theta_s<=pi)){
		if (cimag(myConfig->Layers[N-1].k)<0.0){
			E.theta = 0.0;
			E.phi = 0.0;
		}else{
			TLGF(zN, z_, kN_rho, 'e', 'v', 0, myConfig, &V_e_v, &I);
			TLGF(zN, z_, kN_rho, 'h', 'i', 0, myConfig, &V_h_i, &I);	
			TLGF(zN, z_, kN_rho, 'h', 'v', 0, myConfig, &V_h_v, &I);	
			E_theta =-V_e_v*(-sin(phi_s)*Ilx+cos(phi_s)*Ily);
			E_phi = +V_h_v*cos(theta_s)*(cos(phi_s)*Ilx+sin(phi_s)*Ily)
					-V_h_i*(1.0/etaN)*(muN/mu_)*sin(theta_s)*cos(theta_s)*(Ilz);
			E.theta = E_theta*ExpN;
			E.phi = E_phi*ExpN;
		}
		flag++;
	}
	assert(flag==1);
	return E;
}

complex double GEJ0xx(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta, complex double eps, double omega){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return j/(3.0*omega*eps0*eps);
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*(g+(g1/R+((x-x_)*(x-x_)/(R*R))*(g2-g1/R))/(k*k));
	}
}

complex double GEJ0yy(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta, complex double eps, double omega){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return j/(3.0*omega*eps0*eps);
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*(g+(g1/R+((y-y_)*(y-y_)/(R*R))*(g2-g1/R))/(k*k));
	}
}

complex double GEJ0zz(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta, complex double eps, double omega){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return j/(3.0*omega*eps0*eps);
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*(g+(g1/R+((z-z_)*(z-z_)/(R*R))*(g2-g1/R))/(k*k));
	}
}

complex double GEJ0xy(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*((x-x_)*(y-y_)/(R*R)*(g2-g1/R)/(k*k));
	}
}

complex double GEJ0xz(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*((x-x_)*(z-z_)/(R*R)*(g2-g1/R)/(k*k));
	}
}

complex double GEJ0yx(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*((y-y_)*(x-x_)/(R*R)*(g2-g1/R)/(k*k));
	}
}

complex double GEJ0yz(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*((y-y_)*(z-z_)/(R*R)*(g2-g1/R)/(k*k));
	}
}

complex double GEJ0zx(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*((z-z_)*(x-x_)/(R*R)*(g2-g1/R)/(k*k));
	}
}

complex double GEJ0zy(double x, double x_, double y, double y_, double z, double z_,
	complex double k, complex double eta){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		complex double g2=-(1.0+j*k*R)*g1/R+g/(R*R);
		return -j*k*eta*((z-z_)*(y-y_)/(R*R)*(g2-g1/R)/(k*k));
	}
}

complex double GEM0xx(){
	return 0.0;
}

complex double GEM0xy(double x, double x_, double y, double y_, double z, double z_,
	complex double k){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		return +((z-z_)/R*g1);
	}
}

complex double GEM0xz(double x, double x_, double y, double y_, double z, double z_,
	complex double k){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		return -((y-y_)/R*g1);
	}
}

complex double GEM0yx(double x, double x_, double y, double y_, double z, double z_,
	complex double k){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		return -((z-z_)/R*g1);
	}
}

complex double GEM0yy(){
	return 0.0;
}

complex double GEM0yz(double x, double x_, double y, double y_, double z, double z_,
	complex double k){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		return +((x-x_)/R*g1);
	}
}

complex double GEM0zx(double x, double x_, double y, double y_, double z, double z_,
	complex double k){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		return +((y-y_)/R*g1);
	}
}

complex double GEM0zy(double x, double x_, double y, double y_, double z, double z_,
	complex double k){
	complex double j=csqrt(-1.0);
	double R=sqrt((x-x_)*(x-x_)+
				  (y-y_)*(y-y_)+
				  (z-z_)*(z-z_));
	if (R==0.0){
		return 0.0;
	}else{
		complex double g=cexp(-j*k*R)/(4.0*pi*R);
		complex double g1=-(1.0+j*k*R)*g/R;
		return -((x-x_)/R*g1);
	}
}

complex double GEM0zz(){
	return 0.0;
}

typedef struct Args Args;
struct Args{
	double z, z_, rho;
	double a, b;
	MLayers *myConfig;
};

complex double GTest_integrand(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho=t, d_k_rho=1.0;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((1.0*V_e_i+1.0*V_h_i)*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJxx_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((V_e_i+V_h_i)*besselj(0, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJxx_integrand2(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((V_e_i-V_h_i)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJxx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1, G2;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJxx_integrand1, &myArgs, a, tol, kmax);
	G2 = SI(GEJxx_integrand2, &myArgs, a, tol, kmax);
	if (m==n){
		return -0.5*G1+0.5*cos(2.0*phi)*G2
			+GEJ0xx(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta, myConfig->Layers[m-1].eps, 2.0*pi*myConfig->freq);
	}else{
		return -0.5*G1+0.5*cos(2.0*phi)*G2;
	}
}

complex double GEJxy_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((V_e_i-V_h_i)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJxy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJxy_integrand1, &myArgs, a, tol, kmax);
	if (m==n){
		return +0.5*sin(2.0*phi)*G1
			+GEJ0xy(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta);
	}else{
		return +0.5*sin(2.0*phi)*G1;
	}
}

complex double GEJyx_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((V_e_i-V_h_i)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJyx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJyx_integrand1, &myArgs, a, tol, kmax);
	if (m==n){
		return +0.5*sin(2.0*phi)*G1
			+GEJ0yx(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta);
	}else{
		return +0.5*sin(2.0*phi)*G1;
	}
}

complex double GEJyy_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((V_e_i+V_h_i)*besselj(0, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJyy_integrand2(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_i, V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V_e_i, &I);
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return ((V_e_i-V_h_i)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJyy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1, G2;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJyy_integrand1, &myArgs, a, tol, kmax);
	G2 = SI(GEJyy_integrand2, &myArgs, a, tol, kmax);
	if (m==n){
		return -0.5*G1-0.5*cos(2.0*phi)*G2
			+GEJ0yy(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta, myConfig->Layers[m-1].eps, 2.0*pi*myConfig->freq);
	}else{
		return -0.5*G1-0.5*cos(2.0*phi)*G2;
	}
}

complex double GEJxz_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	return (k_rho*V_e_v*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJxz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJxz_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		return (eta0/(j*k0*eps_))*cos(phi)*G1
			+GEJ0xz(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta);
	}else{
		return (eta0/(j*k0*eps_))*cos(phi)*G1;
	}
}

complex double GEJyz_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	return (k_rho*V_e_v*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJyz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJyz_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps_=myConfig->Layers[m-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		return (eta0/(j*k0*eps_))*sin(phi)*G1
			+GEJ0yz(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta);
	}else{
		return (eta0/(j*k0*eps_))*sin(phi)*G1;
	}
}

complex double GEJzx_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V;
	complex double I_e_i;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V, &I_e_i);
	return (k_rho*I_e_i*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJzx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJzx_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		return (eta0/(j*k0*eps))*cos(phi)*G1
			+GEJ0zx(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta);
	}else{
		return (eta0/(j*k0*eps))*cos(phi)*G1;
	}
}

complex double GEJzy_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V;
	complex double I_e_i;
	TLGF(z, z_, k_rho, 'e', 'i', 1, myConfig, &V, &I_e_i);
	return (k_rho*I_e_i*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJzy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJzy_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		return (eta0/(j*k0*eps))*sin(phi)*G1
			+GEJ0zy(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
					myConfig->Layers[m-1].eta);
	}else{
		return (eta0/(j*k0*eps))*sin(phi)*G1;
	}
}

complex double GEJzz_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V;
	complex double I_e_v;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V, &I_e_v);
	return (k_rho*k_rho*I_e_v*besselj(0, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEJzz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEJzz_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n-1].eps;
	complex double eps_=myConfig->Layers[m-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		if ((rho==0.0)&&(z==z_)){
			return -(eta0*eta0/(k0*k0*eps*eps_))*G1+j*(eta0/(k0*eps))
					+GEJ0zz(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
						myConfig->Layers[m-1].eta, myConfig->Layers[m-1].eps, 2.0*pi*myConfig->freq);
		}else{
			return -(eta0*eta0/(k0*k0*eps*eps_))*G1
					+GEJ0zz(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k, 
						myConfig->Layers[m-1].eta, myConfig->Layers[m-1].eps, 2.0*pi*myConfig->freq);
		}
		
	}else{
		return -(eta0*eta0/(k0*k0*eps*eps_))*G1;
	}
}

complex double GEMxx_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v, V_h_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	TLGF(z, z_, k_rho, 'h', 'v', 1, myConfig, &V_h_v, &I);
	return ((V_e_v-V_h_v)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMxx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMxx_integrand1, &myArgs, a, tol, kmax);
	if (m==n){
		return -0.5*sin(2.0*phi)*G1
			+GEM0xx();
	}else{
		return -0.5*sin(2.0*phi)*G1;
	}
}

complex double GEMxy_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v, V_h_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	TLGF(z, z_, k_rho, 'h', 'v', 1, myConfig, &V_h_v, &I);
	return ((V_e_v+V_h_v)*besselj(0, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMxy_integrand2(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v, V_h_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	TLGF(z, z_, k_rho, 'h', 'v', 1, myConfig, &V_h_v, &I);
	return ((V_e_v-V_h_v)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMxy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1, G2;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMxy_integrand1, &myArgs, a, tol, kmax);
	G2 = SI(GEMxy_integrand2, &myArgs, a, tol, kmax);
	if (m==n){
		return -0.5*G1+0.5*cos(2.0*phi)*G2
			+GEM0xy(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k);
	}else{
		return -0.5*G1+0.5*cos(2.0*phi)*G2;
	}
}

complex double GEMyx_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v, V_h_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	TLGF(z, z_, k_rho, 'h', 'v', 1, myConfig, &V_h_v, &I);
	return ((V_e_v+V_h_v)*besselj(0, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMyx_integrand2(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v, V_h_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	TLGF(z, z_, k_rho, 'h', 'v', 1, myConfig, &V_h_v, &I);
	return ((V_e_v-V_h_v)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMyx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1, G2;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMyx_integrand1, &myArgs, a, tol, kmax);
	G2 = SI(GEMyx_integrand2, &myArgs, a, tol, kmax);
	if (m==n){
		return +0.5*G1+0.5*cos(2.0*phi)*G2
			+GEM0yx(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k);
	}else{
		return +0.5*G1+0.5*cos(2.0*phi)*G2;
	}
}

complex double GEMyy_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_e_v, V_h_v;
	complex double I;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V_e_v, &I);
	TLGF(z, z_, k_rho, 'h', 'v', 1, myConfig, &V_h_v, &I);
	return ((V_e_v-V_h_v)*besselj(2, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMyy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMyy_integrand1, &myArgs, a, tol, kmax);
	if (m==n){
		return +0.5*sin(2.0*phi)*G1
			+GEM0yy();
	}else{
		return +0.5*sin(2.0*phi)*G1;
	}
}

complex double GEMxz_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return (k_rho*V_h_i*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMxz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMxz_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m-1].mu;
	double k0=myConfig->k0;
	if (m==n){
		return +(1.0/(j*eta0*k0*mu_))*sin(phi)*G1
			+GEM0xz(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k);
	}else{
		return +(1.0/(j*eta0*k0*mu_))*sin(phi)*G1;
	}
}

complex double GEMyz_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V_h_i;
	complex double I;
	TLGF(z, z_, k_rho, 'h', 'i', 1, myConfig, &V_h_i, &I);
	return (k_rho*V_h_i*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMyz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMyz_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double mu_=myConfig->Layers[m-1].mu;
	double k0=myConfig->k0;
	if (m==n){
		return -(1.0/(j*eta0*k0*mu_))*cos(phi)*G1
			+GEM0yz(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k);
	}else{
		return -(1.0/(j*eta0*k0*mu_))*cos(phi)*G1;
	}
}

complex double GEMzx_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V;
	complex double I_e_v;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V, &I_e_v);
	return (k_rho*I_e_v*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMzx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMzx_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		return +((j*eta0)/(k0*eps))*sin(phi)*G1
			+GEM0zx(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k);
	}else{
		return +((j*eta0)/(k0*eps))*sin(phi)*G1;
	}
}

complex double GEMzy_integrand1(complex double t, void *args){
	Args *myArgs=(Args*)args;
	double z=myArgs->z;
	double z_=myArgs->z_;
	double rho=myArgs->rho;
	double a=myArgs->a;
	double b=myArgs->b;
	MLayers *myConfig=myArgs->myConfig;
	complex double k_rho, d_k_rho;
	detour_k_rho(t, a, b, &k_rho, &d_k_rho);
	complex double V;
	complex double I_e_v;
	TLGF(z, z_, k_rho, 'e', 'v', 1, myConfig, &V, &I_e_v);
	return (k_rho*I_e_v*besselj(1, k_rho*rho)*k_rho*d_k_rho)/(2.0*pi);
}

complex double GEMzy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	double rho=sqrt((x-x_)*(x-x_)+(y-y_)*(y-y_));
	double phi=atan2(y-y_, x-x_);
	complex double G1;
	double a=findMax(myConfig);
	double b=findDetourParameters(rho, z, z_, myConfig);
	Args myArgs={z, z_, rho, a, b, myConfig};
	int m, n;
	m = selectLayer(z_, myConfig);
	n = selectLayer(z, myConfig);
	G1 = SI(GEMzy_integrand1, &myArgs, a, tol, kmax);
	complex double j=csqrt(-1.0);
	complex double eps=myConfig->Layers[n-1].eps;
	double k0=myConfig->k0;
	if (m==n){
		return -((j*eta0)/(k0*eps))*cos(phi)*G1
			+GEM0zy(x, x_, y, y_, z, z_, myConfig->Layers[m-1].k);
	}else{
		return -((j*eta0)/(k0*eps))*cos(phi)*G1;
	}
}

complex double GEMzz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig){
	assert(myConfig->N>0);
	return 0.0*x*y*z*x_*y_*z_*tol*kmax;
}