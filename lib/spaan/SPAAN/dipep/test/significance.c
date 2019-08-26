#include<stdio.h>
#define I 65
#define M 4000

void check(double in[M][I],int nin,int min)
{
	FILE *fp=fopen("test.txt","w");
	int i,j;
	for(i=0;i<min;i++)
	{
		for(j=0;j<=nin;j++)
			fprintf(fp,"%lf ",in[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return;
}

void significance(int argc,char *argv[])
{
	FILE *fp=fopen(argv[1],"r");
	int i,j,k;
	int n;
	int nin;
	int ad,nonad;
	int signif[I];
	int add[M];
	float limit=1.0;
	double in[M][I];
	double avg[2][I];
	
	printf("Enter the number of input parameters: ");
	scanf("%d",&nin);
	printf("Enter the limit: ");
	scanf("%f",&limit);
	for(i=0;i<=nin;i++)
	{
		avg[0][i]=0;
		avg[1][i]=0;
		signif[i]=0;
	}
	
	for(i=0,k=1;k==1;i++)
	{
		for(j=0;j<nin;j++)
			if(fscanf(fp,"%lf",&in[i][j])==EOF)
			{
				k=0;
				break;
			}
		if(fscanf(fp,"%d",&add[i])==EOF)
		{
			k=0;
			break;
		}
	}
	n=i;


	for(i=0,ad=0,nonad=0;i<n;i++)
	{
		if(add[i]==0)
		{
			nonad++;
			for(j=0;j<nin;j++)	
			{
				avg[0][j]=(avg[0][j]*(nonad-1)+in[i][j])/nonad;
			}
		}
	}


	for(i=0;i<n;i++)
	{
		if(add[i]==1)
		{
			for(j=0;j<nin;j++)
			{
				in[i][j]=in[i][j]/avg[0][j];
				//printf("%lf/%lf ",in[i][j],avg[0][j]);
				if(in[i][j]<=limit)
					signif[j]++;
			}
			//printf(" %dth\n",i);
		}
	}

	ad=n-nonad;
	for(i=0;i<nin;i++)
		printf("%.2f\n",100.0*signif[i]/ad);
	printf("\n");
	fclose(fp);
	check(in,nin,n);
	return;
}

main(int argc,char *argv[])
{
	significance(argc,argv);
}

