#include<stdio.h>
#include<stdlib.h>
#define N 12288
#define KK 600
static double T[N][KK]={0.0},T1[N][KK]={0.0},T2[N][KK]={0.0};
int main()
{
	int i,j;
	FILE *fp;
	fp=fopen("J_a.txt","w");

	FILE *fp1;
	fp1=fopen("J_a1.txt","r");

	FILE *fp2;
	fp2=fopen("J_a2.txt","r");

	for(i=0;i<N;i++)
	{
		for(j=0;j<KK;j++)
		{
			fscanf(fp1,"%lf",&T1[i][j]);
			fscanf(fp2,"%lf",&T2[i][j]);

			T[i][j]=(T1[i][j]+T2[i][j]);
			fprintf(fp,"%f ",T[i][j]);
		}
		fprintf(fp,"\n");
	}

	return 0;

}