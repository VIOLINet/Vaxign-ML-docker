/*-- This program converts the input protein sequence file to a standard format, and saves it in a file with extension '.std' --*/

#include<stdio.h>
#include<math.h>
#include"./path.h"
#define M 100000
#define N 10000

//Pritilipi Adhikaar:
//			 Sachu, Bhootpurve chhatra,
//			 Bhartiya Prodyogiki Sansthan, Kanpur
//			 Dinaank: Prathma, Shukl Paksh,
//			 Asharh Mah, Vikram Samwat 2054
//			 Godhuli Vela, Saymkaal

void prnseq(FILE *fout,char *pro,int np,char *name,int nn)
{
	int i;
	for(i=0;i<nn;i++)
		fprintf(fout,"%c",name[i]);
	for(i=0;i<np;i++)
		fprintf(fout,"%c",pro[i]);
	fprintf(fout,"\n");
	return;
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

void readfile(int argc, char *argv[])
{

	FILE *fin,*fout;
	char temp;
	char infile[50], outfile[50];
//	char pathin[50]="/home/sachu/test/";
//	char pathout[50]="/home/sachu/test/";
	char inext[10]="fasta";
	char outext[10]="std";
	char pro[M],name[N];	//array which stores the sequence of a protein
	int n=0;	//no of protein sequences
	int nfinal=0;
	int i=0,j=0,check=0;
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

	while(fscanf(fin,"%c",&temp)!=EOF)
	{
		if(temp!='>')
		{
			if(temp!='\n')
			{
				pro[i]=temp;
				i++;
			}
		}
		else
		{
			if(n!=0)
				prnseq(fout,pro,i,name,j);
			n++;
			i=0;
			j=0;
			name[j++]=temp;
			do{
				fscanf(fin,"%c",&temp);
				name[j++]=temp;
			}while(temp!='\n');
		}
	}
	prnseq(fout,pro,i,name,j);
	printf("No of seq in the input file : %d\n",n);
	fclose(fin);
	fclose(fout);
}

int main(int argc, char *argv[])
{
	readfile(argc,argv);
	return 0;
}
