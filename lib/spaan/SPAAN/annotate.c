/*-- This program converts the input protein sequence file to a standard format, and saves it in a file with extension '.std' --*/

#include<stdio.h>
#include<math.h>
#include"./path.h"
#define M 100000
#define N 1000000

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

void prnseq_pos(FILE *fin,FILE *fout,int pos)
{
	int i;
	char temp;
	fseek(fin,pos,SEEK_SET);
	for(i=0;i<2;i++)
	{
		while(fscanf(fin,"%c",&temp)!=EOF && temp!='\n')
			fprintf(fout,"%c",temp);
		fprintf(fout,"\n");
	}
	return;
}

int prn_out(FILE *fin,FILE *fout,int *pos,int *out_pro,int total)
{
	int i;
	int no=0;
	for(i=0;i<total;i++)
		if(out_pro[i]!=0)
		{
			prnseq_pos(fin,fout,pos[i]);
			no++;
		}
	return no;
}

int line_length(FILE *fin,int pos)
{
	int length=0;
	char temp;
	fseek(fin,pos,SEEK_SET);
	for(fscanf(fin,"%c",&temp);temp!='\n';fscanf(fin,"%c",&temp));
	while(fscanf(fin,"%c",&temp)!=EOF && temp!='\n')
		length++;
	return length;
}

int initpos(FILE *fin,int *pos)
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

void Out_pro(FILE *fin,FILE *fin1,int *op,int *pos,int total)
{
	int i,no=0;
	int gi[2];
	int first=1;
	int max=0;
	char temp[20];
	int score,x,y,xprev=0;

	for(i=0;i<total;i++)
		op[i]=1;


	while(fscanf(fin1,"%s",temp)!=EOF)
	{
		fseek(fin1,2,SEEK_CUR);
		fscanf(fin1,"%d",&x);
		fseek(fin1,1,SEEK_CUR);
		fscanf(fin1,"%d",&y);
		fseek(fin1,19,SEEK_CUR);
		fscanf(fin1,"%d",&score);
		if(op[x-1]==1)
		{
			if(score >= 90)
			{
			/*	fseek(fin,pos[x-1]+4,SEEK_SET);
				fscanf(fin,"%d",&gi[0]);
				printf("%d ",gi[0]);
				fseek(fin,pos[y-1]+4,SEEK_SET);
				fscanf(fin,"%d",&gi[1]);
				printf("%d %d\n",gi[1],score);
			*/
				if(first==1)
				{
					first=0;
					xprev=x;
					op[y-1]=0;
					if(line_length(fin,pos[x-1])>line_length(fin,pos[y-1]))
						max=x;
					else
						max=y;
				}
				else if(x==xprev)
				{
					if(line_length(fin,pos[y-1])>line_length(fin,pos[max-1]))
						max=y;
					op[y-1]=0;
				}
				else
				{
					op[xprev-1]=0;
					op[max-1]=2;
					xprev=x;
					op[y-1]=0;
					if(line_length(fin,pos[x-1])>line_length(fin,pos[y-1]))
						max=x;
					else
						max=y;
				}
			}
		}
	}

	op[xprev-1]=0;
	op[max-1]=2;

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
	//char pathin[50]="/home/sachu/test/";
	//char pathout[50]="/home/sachu/test/";
	char inext[10]="flt";
	char outext[10]="ant";
	char pro[M],name[N];	//array which stores the sequence of a protein
	char annotation[50];
	int n=0;	//no of protein sequences
	int nfinal=0;
	int i=0,j=0,check=0;
	int total;
	int pos[N];
	if(argc!=3)
	{
		printf("There must be two arguments. (Input file name & annotation)\n");
		printf("Enter the inputfile: ");
		scanf("%s",infile);
		printf("Annotation: ");
		scanf("%s",annotation);
	}
	else
	{
		strcpy(infile,argv[1]);
		strcpy(annotation,argv[2]);
	}

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

	total=initpos(fin,pos);

	for(i=0;i<total;i++)
	{
		fseek(fin,pos[i]+1,SEEK_SET);
		fprintf(fout,">%s",annotation);
		while(fscanf(fin,"%c",&temp)!=EOF && temp!='\n')
			fprintf(fout,"%c",temp);
		fprintf(fout,"\n");
		while(fscanf(fin,"%c",&temp)!=EOF && temp!='\n')
				fprintf(fout,"%c",temp);
		fprintf(fout,"\n");
	}

	printf("No of seq in the input file : %d\n",total);
	fclose(fin);
	fclose(fout);
}

int main(int argc, char *argv[])
{
	readfile(argc,argv);
	return 0;
}
