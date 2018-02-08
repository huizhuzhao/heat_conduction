#include<stdio.h>
#include<stdlib.h>
#include"mpi.h"
#include<math.h>
#include<time.h>
#define pi 3.141593
int B_C=2;
#define N 0
#define N1 6144   // Lattice 两端热库加热部分的各自长度
#define N_tot 12288  // Lattice 的总长度

#define KK 600  
#define st 160
#define N_size 64

static double e=0.1,E0=e*N_tot,T0=1.0,A=0.2,Z=1024.0,Tt=47.0*3600.0;  // e是Lattice的能量密度 Tt是程序的运行时间
static int M_chi=100000,M_yu=20000; // M_chi 是Lattice第一次运行的预热时间  M_yu是Lattice在前一次的温度基础上，以后每次的加热时间
static int M_xun[N_size]={0},M_xun1[N_size]={0}; // M_xun[rank]记录编号为 rank 的核运行的循环次数， M_xun1[]将所有核的循环次数通过MPI_REDUCE函数收集起来
static double h=0.01,q[N_tot]={0.0},p[N_tot]={0.0},T_t[N_tot][KK]={0.0},E_e[N_tot][KK]={0.0},J_a[N_tot][KK]={0.0},
T_t1[N_tot][KK]={0.0},E_e1[N_tot][KK]={0.0},J_a1[N_tot][KK]={0.0};

double r=0.0,z[N1*2]={0.0},t[2*N1]={0.0};///Langevin

int main(int argc,char **argv)
{
	int i,j,k,rank,size;
	MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

	FILE *ft;
	ft=fopen("time.txt","w");
	
	FILE *fp_E;
	fp_E=fopen("E.txt","w");

	FILE *fp_T;
	fp_T=fopen("T.txt","w");

	FILE *fp_J_a;
	fp_J_a=fopen("J_a.txt","w");




	time_t start,end;
	double diff=0.0;
	srand(time(NULL)+200*rank);////////////*******************************
	start=time(NULL);

	double v1(double x);
	double v(double x);
	double U(double x);
	double F(double x);
	void RK_4();
	void RK6_Butcher();
	void RKN_89_2();
	double gauss();
	double E_fixed();
	double E_periodic();
	double E_free();
	double Momentum();
	void E_T_J(double *q,double *p,double *E,double *T,double *J_a);  ////E[N_tot]; T[N_tot]; J[N_tot];
    

	
	for(k=0;k<N1;k++)
	{
		t[k]=T0+A*cos(2.0*pi*k/Z);t[k+N1]=T0+A*cos(2.0*pi*(k+N1)/Z);
	}
	double tem_E[N_tot]={0.0},tem_T[N_tot]={0.0},tem_J_a[N_tot]={0.0};

	for(k=0;k<N_tot;k++)
		{
			q[k]=0.0;p[k]=0.0;
		}
	for(i=0;i<M_chi;i++)
	{
		h=0.01; r=0.1;
			for(k=0;k<N1;k++)
			{
				z[k]=gauss()*sqrt(2.0*r*t[k]/h);z[k+N1]=gauss()*sqrt(2.0*r*t[k+N1]/h);
			}
			RK6_Butcher();
			//RKN_89_2();
	}
	double tem_q[N_tot]={0.0},tem_p[N_tot]={0.0};
	for(k=0;k<N_tot;k++)
	{
		tem_q[k]=q[k];tem_p[k]=p[k];
	}
	for(i=0;diff<Tt;i++)
	{
		for(k=0;k<N_tot;k++)
		{
			q[k]=tem_q[k];p[k]=tem_p[k];
		}
		for(j=0;j<M_yu;j++)
		{
			h=0.01; r=0.1;
			for(k=0;k<N1;k++)
			{
				z[k]=gauss()*sqrt(2.0*r*t[k]/h);z[k+N1]=gauss()*sqrt(2.0*r*t[k+N1]/h);
			}
			RK6_Butcher();
			//RKN_89_2();
		}
		for(k=0;k<N_tot;k++)
		{
		    tem_q[k]=q[k];tem_p[k]=p[k];
		}

		int kk=0;
	    kk=0;
	    for(i=0;i<KK*st;i++)
		{
			if(i%st==0)
			{
				E_T_J(q,p,tem_E,tem_T,tem_J_a);

				for(k=0;k<N_tot;k++)
				{
					E_e[k][kk]+=tem_E[k];
					T_t[k][kk]+=tem_T[k];
					J_a[k][kk]+=tem_J_a[k];
					//J_b[k][kk]+=tem_J_b[k];
					//J_c[k][kk]+=tem_J_c[k];
				}
				
				kk++;
			}
			h=0.05; r=0.0;
			for(k=0;k<N1;k++)
			{
				z[k]=0.0;//gauss()*sqrt(2.0*r*t[k]/h);
				z[k+N1]=0.0;//gauss()*sqrt(2.0*r*t[k+N1]/h);
			}
			//RK6_Butcher();
			RKN_89_2();
		}
		    
		M_xun[rank]++;

	end=time(NULL);
	diff=difftime(end,start);
	}
	
	MPI_Reduce(&T_t[0][0],&T_t1[0][0],N_tot*KK,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&E_e[0][0],&E_e1[0][0],N_tot*KK,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&J_a[0][0],&J_a1[0][0],N_tot*KK,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//	MPI_Reduce(&J_b[0][0],&J_b1[0][0],N_tot*KK,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//	MPI_Reduce(&J_c[0][0],&J_c1[0][0],N_tot*KK,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&M_xun[0],&M_xun1[0],N_size,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

	if(rank==0)
	{
		for(k=0;k<N_tot;k++)
		{
			for(j=0;j<KK;j++)
			{
				fprintf(fp_T,"%f ",T_t1[k][j]); fprintf(fp_E,"%f ",E_e1[k][j]);
				fprintf(fp_J_a,"%f ",J_a1[k][j]); //fprintf(fp_J_b,"%f ",J_b1[k][j]);
			}
			fprintf(fp_T,"\n");fprintf(fp_E,"\n");
			fprintf(fp_J_a,"\n"); //fprintf(fp_J_b,"\n"); 
		}

		int tem_tot=0;
		for(i=0;i<size;i++)
		    tem_tot+=M_xun1[i];

		fprintf(ft,"time: %f sec    %f min   %f hour   BC:%d\n\n",diff,diff/60.0,diff/3600.0,B_C);
		fprintf(ft,"L=%d Z=%d T0=%f  A=%f  e=%f  KK=%d st=%d  h=%f\n",N_tot,(int)Z,T0,A,e,KK,st,h);
		fprintf(ft,"size=%d tem_tot=%d\n\n\n",size,tem_tot);

		for(i=0;i<size;i++)
			fprintf(ft,"%d\n",M_xun1[i]);
	}


    

    MPI_Finalize();
    fclose(fp_E);
	fclose(fp_T);
    fclose(fp_J_a);
    fclose(ft);


    return 0;
}
double v(double x)
{
	double y;
	y=x*x/2.0+x*x*x*2.0/3.0+x*x*x*x/4.0;
	return y;
}
double v1(double x)
{
	double y;
	y=x+2.0*x*x+x*x*x;
	return y;
}
double U(double x)//////外势
{
	double y;
	y=0.0;//10.0*x*x*x*x/4.0;
	return y;
}
double F(double x)
{
	double y;
	y=0.0;//10.0*x*x*x;
	return y;
}
void kf(double *x,double *y)
{
	int i;
	for(i=0;i<N_tot;i++)
		y[i]=x[i];
}
void wf(double *q,double *p,double *y) //需要 N_tot,r,z[]
{
	int i;
	double f[N_tot+1]={0.0};
	double v1(double);

	switch (B_C)  //1: Fixed BC  2: Periodic BC  3: Free BC
	{
	      case 1: f[0]=-v1(q[0]);f[N_tot]=-v1(-q[N_tot-1]);  //Fixed BC
			               break;
          case 2: f[0]=-v1(q[0]-q[N_tot-1]);f[N_tot]=f[0];  //Periodic BC
			               break;
          case 3: f[0]=0.0; f[N_tot]=0.0;    //Free BC
			               break;

	}
	for(i=1;i<N_tot;i++)
	{
		f[i]=-v1(q[i]-q[i-1]);
	}

	for(i=0;i<N_tot;i++)
	{
		y[i]=f[i]-f[i+1]-F(q[i]);
	}
	
	for(i=0;i<N1;i++)                                   //// Langevin 热库
	{
		y[i]=y[i]-r*p[i]+z[i];
		y[N_tot-N1+i]=y[N_tot-N1+i]-r*p[N_tot-N1+i]+z[N1+i];
	}


}
void RK_4()  //需要 h,N_tot,q[],p[],kf(),wf() 
{
	int i;
	static double tem_q[N_tot]={0.0},tem_p[N_tot]={0.0};
	static double k1[N_tot]={0.0},k2[N_tot]={0.0},k3[N_tot]={0.0},k4[N_tot]={0.0},
		          w1[N_tot]={0.0},w2[N_tot]={0.0},w3[N_tot]={0.0},w4[N_tot]={0.0};
	void kf(double *x,double *y);
	void wf(double *q,double *p,double *y);

	
    kf(p,k1);
	wf(q,p,w1);
	
	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h/2.0*k1[i];
		tem_p[i]=p[i]+h/2.0*w1[i];
	}
	kf(tem_p,k2);
	wf(tem_q,tem_p,w2);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h/2.0*k2[i];
		tem_p[i]=p[i]+h/2.0*w2[i];
	}
	kf(tem_p,k3);
	wf(tem_q,tem_p,w3);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h*k3[i];
		tem_p[i]=p[i]+h*w3[i];
	}
	kf(tem_p,k4);
	wf(tem_q,tem_p,w4);
    
	for(i=0;i<N_tot;i++)
	{
		q[i]=q[i]+h/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
		p[i]=p[i]+h/6.0*(w1[i]+2.0*w2[i]+2.0*w3[i]+w4[i]);
	}

}
void RK6_Butcher() //需要 h,N_tot,q[],p[],kf(),wf() 
{
	int i;
    static double tem_q[N_tot]={0.0},tem_p[N_tot]={0.0};
	static double k1[N_tot]={0.0},k2[N_tot]={0.0},k3[N_tot]={0.0},k4[N_tot]={0.0},k5[N_tot]={0.0},k6[N_tot]={0.0},
		          w1[N_tot]={0.0},w2[N_tot]={0.0},w3[N_tot]={0.0},w4[N_tot]={0.0},w5[N_tot]={0.0},w6[N_tot]={0.0};
	void kf(double *x,double *y);
	void wf(double *q,double *p,double *y);

	
    kf(p,k1);
	wf(q,p,w1);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h/4.0*k1[i];
		tem_p[i]=p[i]+h/4.0*w1[i];
	}
	kf(tem_p,k2);
	wf(tem_q,tem_p,w2);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h/8.0*(k1[i]+k2[i]);
		tem_p[i]=p[i]+h/8.0*(w1[i]+w2[i]);
	}
	kf(tem_p,k3);
	wf(tem_q,tem_p,w3);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h*(-k2[i]/2.0+k3[i]);
		tem_p[i]=p[i]+h*(-w2[i]/2.0+w3[i]);
	}
	kf(tem_p,k4);
	wf(tem_q,tem_p,w4);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h*(3.0*k1[i]+9.0*k4[i])/16.0;
		tem_p[i]=p[i]+h*(3.0*w1[i]+9.0*w4[i])/16.0;
	}
	kf(tem_p,k5);
	wf(tem_q,tem_p,w5);

	for(i=0;i<N_tot;i++)
	{
		tem_q[i]=q[i]+h*(-3.0*k1[i]+2.0*k2[i]+12.0*k3[i]-12.0*k4[i]+8.0*k5[i])/7.0;
		tem_p[i]=p[i]+h*(-3.0*w1[i]+2.0*w2[i]+12.0*w3[i]-12.0*w4[i]+8.0*w5[i])/7.0;
	}
	kf(tem_p,k6);
	wf(tem_q,tem_p,w6);

	for(i=0;i<N_tot;i++)
	{
		q[i]=q[i]+h/90.0*(7.0*k1[i]+32.0*k3[i]+12.0*k4[i]+32.0*k5[i]+7.0*k6[i]);
		p[i]=p[i]+h/90.0*(7.0*w1[i]+32.0*w3[i]+12.0*w4[i]+32.0*w5[i]+7.0*w6[i]);
	}

}

void RKN_89() //需要 h,N_tot,q[],p[],kf(q,y),wf(q,p,y) 
{
   double c[9]={0.0, 1.0/20.0, 1.0/10.0, 3.0/10.0, 1.0/2.0, 7.0/10.0, 9.0/10.0, 1.0, 1.0},
	   b[9]={7987313.0/109941300.0, 0.0, 1610737.0/44674560.0, 10023263.0/33505920.0, -497221.0/12409600.0, 10023263.0/78180480.0, 1610737.0/402071040.0, 0.0, 0.0},
       b1[9]={7987313.0/109941300.0, 0.0, 1610737.0/40207104.0, 10023263.0/23454144.0, -497221.0/6204800.0, 10023263.0/23454144.0, 1610737.0/40207104.0, -4251941.0/54970650.0, 3.0/20.0},
	   b_1[9]={223.0/7938.0, 0.0, 1175.0/8064.0, 925.0/6048.0, 41.0/448.0, 925.0/14112.0, 1175.0/72576.0, 0.0, 0.0},
	   b_1_1[9]={223.0/7938.0, 0.0, 5875.0/36288.0, 4625.0/21168.0, 41.0/224.0, 4625.0/21168.0, 5875.0/36288.0, 223.0/7938.0, 0.0};

   double a[9][9]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                1.0/800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				1.0/600.0, 1.0/300.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				9.0/200.0, -9.0/100.0, 9.0/100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				-66701.0/197352.0, 28325.0/32892.0, -2665.0/5482.0, 2170.0/24669.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				227015747.0/304251000.0, -54897451.0/30425100.0, 12942349.0/10141700.0, -9499.0/304251.0, 539.0/9250.0, 0.0, 0.0, 0.0, 0.0,
				-1131891597.0/901789000.0, 41964921.0/12882700.0, -6663147.0/3220675.0, 270954.0/644135.0, -108.0/5875.0, 114.0/1645.0, 0.0, 0.0, 0.0,
				13836959.0/3667458.0, -17731450.0/1833729.0, 1063919505.0/156478208.0, -33213845.0/39119552.0, 13335.0/28544.0, -705.0/14272.0, 1645.0/57088.0, 0.0, 0.0,
				223.0/7938.0, 0.0, 1175.0/8064.0, 925.0/6048.0, 41.0/448.0, 925.0/14112.0, 1175.0/72576.0, 0.0, 0.0 };
   
   
   double v1(double x);
   double F(double x);
   void wf(double *q,double *p,double *y);

   double tem_q[N_tot]={0.0},tem_p[N_tot]={0.0},f1[N_tot+1]={0.0},f[9][N_tot]={0.0};
   int i,j,k,s=9;
   for(i=0;i<s;i++)
   {
	   for(k=0;k<N_tot;k++)
	   {
		   tem_q[k]=q[k]+p[k]*c[i]*h;
		   
		   for(j=0;j<i;j++)
			   tem_q[k]=tem_q[k]+h*h*a[i][j]*f[j][k];//q

	   }
	   wf(tem_q,tem_p,f[i]);
	  
   }
   for(k=0;k<N_tot;k++)
	{
	   q[k]=q[k]+h*p[k];
		for(i=0;i<s;i++)
		{
			q[k]=q[k]+h*h*b[i]*f[i][k];
			p[k]=p[k]+h*b1[i]*f[i][k];
		}
	}
}
void RKN_89_2() //需要 h,N_tot,q[],p[],kf(q,y),wf(q,p,y) 精度更高
{
   double c[9]={0.0, 1.0/20.0, 1.0/10.0, 3.0/10.0, 1.0/2.0, 7.0/10.0, 9.0/10.0, 1.0, 1.0},
	   b[9]={7987313.0/109941300.0, 0.0, 1610737.0/44674560.0, 10023263.0/33505920.0, -497221.0/12409600.0, 10023263.0/78180480.0, 1610737.0/402071040.0, 0.0, 0.0},
       b1[9]={7987313.0/109941300.0, 0.0, 1610737.0/40207104.0, 10023263.0/23454144.0, -497221.0/6204800.0, 10023263.0/23454144.0, 1610737.0/40207104.0, -4251941.0/54970650.0, 3.0/20.0},
	   b_1[9]={223.0/7938.0, 0.0, 1175.0/8064.0, 925.0/6048.0, 41.0/448.0, 925.0/14112.0, 1175.0/72576.0, 0.0, 0.0},
	   b_1_1[9]={223.0/7938.0, 0.0, 5875.0/36288.0, 4625.0/21168.0, 41.0/224.0, 4625.0/21168.0, 5875.0/36288.0, 223.0/7938.0, 0.0};

   double a[9][9]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                1.0/800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				1.0/600.0, 1.0/300.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				9.0/200.0, -9.0/100.0, 9.0/100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				-66701.0/197352.0, 28325.0/32892.0, -2665.0/5482.0, 2170.0/24669.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				227015747.0/304251000.0, -54897451.0/30425100.0, 12942349.0/10141700.0, -9499.0/304251.0, 539.0/9250.0, 0.0, 0.0, 0.0, 0.0,
				-1131891597.0/901789000.0, 41964921.0/12882700.0, -6663147.0/3220675.0, 270954.0/644135.0, -108.0/5875.0, 114.0/1645.0, 0.0, 0.0, 0.0,
				13836959.0/3667458.0, -17731450.0/1833729.0, 1063919505.0/156478208.0, -33213845.0/39119552.0, 13335.0/28544.0, -705.0/14272.0, 1645.0/57088.0, 0.0, 0.0,
				223.0/7938.0, 0.0, 1175.0/8064.0, 925.0/6048.0, 41.0/448.0, 925.0/14112.0, 1175.0/72576.0, 0.0, 0.0 };
   
   
   double v1(double x);
   double F(double x);
   void wf(double *q,double *p,double *y);

   static double tem_q[N_tot]={0.0},tem_p[N_tot]={0.0},f1[N_tot+1]={0.0},f[9][N_tot]={0.0};
   int i,j,k,s=9;
   for(i=0;i<s;i++)
   {
	   for(k=0;k<N_tot;k++)
	   {
		   tem_q[k]=q[k]+p[k]*c[i]*h;
		   
		   for(j=0;j<i;j++)
			   tem_q[k]=tem_q[k]+h*h*a[i][j]*f[j][k];//q

	   }
	   wf(tem_q,tem_p,f[i]);
	  
   }
   for(k=0;k<N_tot;k++)
	{
	   q[k]=q[k]+h*p[k];
		for(i=0;i<s;i++)
		{
			q[k]=q[k]+h*h*b_1[i]*f[i][k];   ///b_1[]   b_1_1[]  阶数更高的
			p[k]=p[k]+h*b_1_1[i]*f[i][k];
		}
	}
}
double E_periodic()
{
	double v(double x);
	double U(double x);
	double y=0.0;
	int i;
	for(i=0;i<N_tot;i++)
		y=y+p[i]*p[i]/2.0;
	for(i=0;i<N_tot-1;i++)
		y=y+v(q[i+1]-q[i]);
	y=y+v(q[0]-q[N_tot-1]);

	for(i=0;i<N_tot;i++)
		y=y+U(q[i]);            //+on site potential
	return y;
}
double E_fixed()
{
	double v(double x);
	double U(double x);
	double y=0.0;
	int i;
	for(i=0;i<N_tot;i++)
		y=y+p[i]*p[i]/2.0;
	for(i=0;i<N_tot-1;i++)
		y=y+v(q[i+1]-q[i]);
	y=y+v(q[0])+v(-q[N_tot-1]);

	for(i=0;i<N_tot;i++)
		y=y+U(q[i]);            //+on site potential
	return y;
}
double E_free()
{
	double v(double x);
	double U(double x);
	double y=0.0;
	int i;
	for(i=0;i<N_tot;i++)
		y=y+p[i]*p[i]/2.0;
	for(i=0;i<N_tot-1;i++)
		y=y+v(q[i+1]-q[i]);

	for(i=0;i<N_tot;i++)
		y=y+U(q[i]);            //+on site potential
	return y;
}
double Momentum()
{
	double y=0.0;
	int i;
	for(i=0;i<N_tot;i++)
		y=y+p[i];
	return y;
}
double gauss()
{
	int j;
	double u1=0.0,u2=0.0,s=0.0,z;
	for(j=1;j<2;j++)
	{
		u1=(rand()%2000*0.001-1);
		u2=(rand()%2000*0.001-1);
		s=u1*u1+u2*u2;
		if(s>0.0&&s<1.0)
		{
		   z=u1*sqrt(-2*log(s)/s);                
		}
		else
			 j--;
	}
	return z;
}

//////////////////
			void E_T_J(double *q,double *p,double *E,double *T,double *J_a)  ////E[N_tot]; T[N_tot]; J[N_tot];
			{
				
				double v(double x);double v1(double x);
				double U(double x);double F(double x);
				int k;

                for(k=0;k<N_tot;k++)
				T[k]=p[k]*p[k];

				for(k=1;k<N_tot-1;k++)
				{
				    E[k]=p[k]*p[k]/2.0+(v(q[k]-q[k-1])+v(q[k+1]-q[k]))/2.0+U(q[k]);
				   // J_a[k]=-p[k]*(v1(q[k-1]-q[k])+v1(q[k]-q[k+1]))/2.0;
					J_a[k]=-(p[k+1]+p[k])*v1(q[k+1]-q[k])/2.0;
			    	//J_b[k]=-p[k]*(v1(q[k]-q[k-1])+v1(q[k+1]-q[k]))/2.0;
				    
				}

				k=0;
					switch (B_C)
				{
			         case 1:
				     E[k]=p[k]*p[k]/2.0+(v(q[k]-0.0)+v(q[k+1]-q[k]))/2.0+U(q[k]);  //Fixed
				    // J_a[k]=-p[k]*(v1(0.0-q[k])+v1(q[k]-q[k+1]))/2.0;
					 J_a[k]=-(p[k+1]+p[k])*v1(q[k+1]-q[k])/2.0;
				     //J_b[k]=-p[k]*(v1(q[k]-0.0)+v1(q[k+1]-q[k]))/2.0;
				     
				     break;

			         case 2:
				     E[k]=p[k]*p[k]/2.0+(v(q[k]-q[N_tot-1])+v(q[k+1]-q[k]))/2.0+U(q[k]);  //Periodic
				    // J_a[k]=-p[k]*(v1(q[N_tot-1]-q[k])+v1(q[k]-q[k+1]))/2.0;
					 J_a[k]=-(p[k+1]+p[k])*v1(q[k+1]-q[k])/2.0;
				     //J_b[k]=-p[k]*(v1(q[k]-q[N_tot-1])+v1(q[k+1]-q[k]))/2.0;
				     
				     break;

			         case 3:
				     E[k]=p[k]*p[k]/2.0+(0.0+v(q[k+1]-q[k]))/2.0+U(q[k]);  //Free
				    // J_a[k]=-p[k]*(0.0+v1(q[k]-q[k+1]))/2.0;
					 J_a[k]=-(p[k+1]+p[k])*v1(q[k+1]-q[k])/2.0;
				     //J_b[k]=-p[k]*(0.0+v1(q[k+1]-q[k]))/2.0;
				     
				     break;
				}
					k=N_tot-1;
				switch (B_C)
				{
			    case 1:
				     E[k]=p[k]*p[k]/2.0+(v(q[k]-q[k-1])+v(0-q[k]))/2.0+U(q[k]);   //Fixed
				    // J_a[k]=-p[k]*(v1(q[k-1]-q[k])+v1(q[k]-0))/2.0;
					 J_a[k]=-(0.0+p[k])*v1(0.0-q[k])/2.0;
				     //J_b[k]=-p[k]*(v1(q[k]-q[k-1])+v1(0-q[k]))/2.0;
				     
				     break;
				
			    case 2: 
				     E[k]=p[k]*p[k]/2.0+(v(q[k]-q[k-1])+v(q[0]-q[k]))/2.0+U(q[k]);  //Periodic
				    // J_a[k]=-p[k]*(v1(q[k-1]-q[k])+v1(q[k]-q[0]))/2.0;
					 J_a[k]=-(p[0]+p[k])*v1(q[0]-q[k])/2.0;
				     //J_b[k]=-p[k]*(v1(q[k]-q[k-1])+v1(q[0]-q[k]))/2.0;
				     
				     break;

			    case 3:
				     E[k]=p[k]*p[k]/2.0+(v(q[k]-q[k-1])+0.0)/2.0+U(q[k]);   //Free
				    // J_a[k]=-p[k]*(v1(q[k-1]-q[k])+0.0)/2.0;
					 J_a[k]=0.0;
				     //J_b[k]=-p[k]*(v1(q[k]-q[k-1])+0.0)/2.0;
				     
                     break;
				}

			}

			//////////////////

		
