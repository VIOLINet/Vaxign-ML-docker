/*--- This program filters out the proteins containing ambiguous amino acids and proteins of length less than 50 AAs. ----*/
#include<stdio.h>
#include<math.h>
#include"./path.h"
#define M 100000
#define N 10000
FILE *fin,*fout;
//Pritilipi Adhikaar:
//			 Sachu, Bhootpurve chhatra,
//			 Bhartiya Prodyogiki Sansthan, Kanpur
//			 Dinaank: Prathma, Shukl Paksh,
//			 Asharh Mah, Vikram Samwat 2054
//			 Godhuli Vela, Saymkaal


int filter(char *pro,int np,char *name,int nn,int nfinal)
{	
	int prn=1;
	int i=0;
	
	for(i=0;i<np;i++)
	{
		if(!(pro[i]=='R'||pro[i]=='K'||pro[i]=='E'||pro[i]=='D'||pro[i]=='S'||pro[i]=='T'||pro[i]=='N'||pro[i]=='Q'||pro[i]=='P'||pro[i]=='H'||pro[i]=='A'||pro[i]=='G'||pro[i]=='Y'||pro[i]=='C'||pro[i]=='W'||pro[i]=='L'||pro[i]=='V'||pro[i]=='I'||pro[i]=='F'||pro[i]=='M'))
			prn=0;
	}
	if(np<50)
		prn=-1;

	if(prn==1)
	{
		nfinal++;
		fprintf(fout,">");
		for(i=0;i<nn;i++)
			fprintf(fout,"%c",name[i]);
		for(i=0;i<np;i++)
			fprintf(fout,"%c",pro[i]);
		fprintf(fout,"\n");
	}

	return nfinal;
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

	char temp;
	char infile[50], outfile[50];
	//char pathin[50]="/home/sachu/test/";
	//char pathout[50]="/home/sachu/test/";
	char inext[10]="ext";
	char outext[10]="flt";
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
				nfinal=filter(pro,i,name,j,nfinal);
			n++;
			i=0;
			j=0;
			do{
				fscanf(fin,"%c",&temp);
				name[j++]=temp;
			}while(temp!='\n');
		}
	}
	nfinal=filter(pro,i,name,j,nfinal);
	printf("No of seq in the original file : %d\n",n);
	printf("No of seq in the filtered file : %d\n",nfinal);
	fclose(fin);
	fclose(fout);
}

int main(int argc, char *argv[])
{
	readfile(argc,argv);
	return 0;
}
