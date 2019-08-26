#include<stdio.h>
#include<time.h>
#include"./path.h"
#include<stdlib.h>
#define M 1000000
#define N 1000000

struct protein{
	char nm[N];
	char sq[M];
	int lnm,lsq;
};
typedef struct protein protein;
//-----Write no functions above it---------

int initpos(int *pos,FILE *fin)
{
	int i;
	int pin;
	char temp;
	int total=0;
	//fseek(fin,0,SEEK_SET);
	while(fscanf(fin,"%c",&temp)!=EOF)
	{
		if(temp=='>')
		{
			pos[total]=ftell(fin)-1;
			total++;
			for(i=0;i<2;i++)
				while(fscanf(fin,"%c",&temp)!=EOF && temp!='\n');
		}
	}

	return total;
}

char *extension(char *a,char *inext,char *outext)
{
	int i;
	for(i=strlen(a)-1;i>=0;i--)
	{
		if(a[i]=='.')
		{
			if(strcmp(&a[i+1],inext)!=0)
			{
				printf("Error: file must have extension .%s\n",inext);
				exit(1);
			}
			a[i+1]='\0';
			break;
		}
	}
	if(i==-1)
	{

		printf("Error: file must have extension .%s\n",inext);
		exit(2);
	}
	strcat(a,outext);
	return a;
}

void prnoutfile(int *pos,int n,FILE *fin,FILE *fout)
{
	int i,j;
	char temp;

	for(i=0;i<n;i++)
	{
		fseek(fin,pos[i],SEEK_SET);
		for(j=0;j<2;j++)
		{
			while(fscanf(fin,"%c",&temp)!=EOF && temp!='\n')
				fprintf(fout,"%c",temp);
			fprintf(fout,"%c",temp);
		}
	}
	return;

}

void shuffle(int *pos,int total)
{
	int i,j,s,temp;
	srand(time(NULL));
	for(i=1;i<=total;i++)
	{
		s=i+fabs(0.5 + ((total-i)*(double)rand()/RAND_MAX));
		temp=pos[i-1];
		pos[i-1]=pos[s-1];
		pos[s-1]=temp;
	}
	return;
}


main(int argc, char *argv[])
{
	FILE *fin, *fout;
	char temp;
	char infile[50], outfile[50];
	//char pathin[50]="/home/sachu/test/";
	//char pathout[50]="/home/sachu/test/";
	char inext[10]="tmp";
	char outext[10]="shf";
	char pro[M],name[N];	//array which stores the sequence of a protein
	int n=0;	//no of protein sequences
	int nfinal=0;
	int i=0,j=0,check=0;
	int pos[N];
	int total=0;
	int nout=0;

	if(argc!=2)
	{
		printf("Enter the inputfile: ");
		scanf("%s",infile);
	}
	else
		strcpy(infile,argv[1]);

	strcpy(outfile,infile);
	extension(outfile,inext,outext);

	strcat(pathin,infile);
	strcat(pathout,outfile);
	fin=fopen(pathin,"r");
	if(fin==NULL)
	{
		printf("File %s not found.\n",infile);
		exit(3);
	}
	fout=fopen(pathout,"w");

	total=initpos(pos,fin);
	printf("Total no of proteins are: %d\n",total);
	shuffle(pos,total);
	printf("Enter the number of proteins to be printed: ");
	scanf("%d",&nout);
	//nout=total;
	prnoutfile(pos,nout,fin,fout);
}
