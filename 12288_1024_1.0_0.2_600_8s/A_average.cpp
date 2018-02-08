#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi 3.141593
#define N 12288
#define Z 1024
#define KK 600
static double T1[N][KK]={0.0},T[Z][KK]={0.0},A[KK]={0.0};

int main()
{
	int i,j,k,m=N/Z;

    FILE *fp1;
	fp1=fopen("12288_1024_1.0_0.2_600_8s_fpu_2ab_J_a.txt","r");  ///原始文件

	FILE *fp;
	fp=fopen("1024_12288_1.0_0.2_600_8s_fpu_2ab_J_a.txt","w");  ///平均后的文件

	FILE *fp_A;
	fp_A=fopen("A_1024_12288_1.0_0.2_600_8s_fpu_2ab_J_a.txt","w");  ///振幅文件


	for(i=0;i<N;i++)
	{
		for(j=0;j<KK;j++)
			fscanf(fp1,"%lf",&T1[i][j]);

	}

	for(i=0;i<Z;i++)
	{
		for(j=0;j<KK;j++)
		{
			for(k=0;k<m;k++)
			{
				T[i][j]=T[i][j]+T1[i+k*Z][j];

			}
			T[i][j]=T[i][j]/(double)m;

			fprintf(fp,"%f ",T[i][j]);
		}
		fprintf(fp,"\n");
	}
	

	double a=0.0,b=0.0;
	double cosi[N]={0.0};
	for(k=0;k<Z;k++)
	{
		cosi[k]=sin(2.0*pi*k/(double)Z);
		b=b+cosi[k]*cosi[k];
	}

	
	for(i=0;i<KK;i++)
	{
		a=0.0;
		for(k=0;k<Z;k++)
			a=a+(T[k][i])*cosi[k];

		A[i]=a/b;

		fprintf(fp_A,"%f ",A[i]);

	}
	

	

	return 0;
}
