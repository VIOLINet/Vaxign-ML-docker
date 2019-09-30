#include<stdio.h>
#include<math.h>
#define M 15000
#define N 15000
#define D 5
#define wmax 75
#define wmin 20
#define gap 10
#define MAX_freq 15000
int max_freq=0;

//Pritilipi Adhikaar:
//			 Sachu, Bhootpurve chhatra,
//			 Bhartiya Prodyogiki Sansthan, Kanpur
//			 Dinaank: Prathma, Shukl Paksh,
//			 Asharh Mah, Vikram Samwat 2054
//			 Godhuli Vela, Saymkaal

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

void freq(int *pro1,int np,int na,char *name)
{
	int i;
	double f[20];
	for(i=0;i<na;i++)
		f[i]=0;
	for(i=0;i<np;i++)
		f[pro1[i]-1]++;
	for(i=0;i<na;i++)
	{
		if(f[i]>max_freq)
			max_freq=f[i];
		//fprintf(fout,"%lf ",f[i]/MAX_freq);
	}
	return;
}
void stats(float *s,int ns,double *m, int order)
{
	int i,j;
	for(i=0;i<order;i++)
		m[i]=0;
	if(ns!=0)
	{	
		for(i=0;i<ns;i++)
			m[0]+=s[i];
		m[0]=m[0]/ns;
		for(j=1;j<order;j++)	//Its starting from 2nd order moment.
		{
			for(i=0;i<ns;i++)
				m[j]+=pow(s[i]-m[0],j+1);
			m[j]=m[j]/ns;
		}
	}
	return;
}

void procalc(FILE *fout,char *pro,int np,char *name,int nn)
{
	int na=5;
	int i=0,j=0,y=0;
	int pro1[M];
	int c=0,k=0,cw=0,cprev=0,check=0,w,wprev;
	int f[5];
	int order=4;
	int clu_pos[100];
	int clu_len[100];
	float pos[5][MAX_freq];
	float pos1[MAX_freq];
	float clu_t[100];
	float f1=0,t,T,tprev;
	double f_clu=0;
	double m[4];
	for(i=0;i<5;i++)
		f[i]=0;

	for(i=0;i<np;i++)
	{
		if(pro[i]=='R'||pro[i]=='K'||pro[i]=='E'||pro[i]=='D')
		{
			pro1[i]=1;
			pos[0][f[0]]=(float)(i+1)/np;
			f[0]++;
		}
		else if(pro[i]=='S'||pro[i]=='T'||pro[i]=='N'||pro[i]=='Q')
		{
			pro1[i]=2;
			pos[1][f[1]]=(float)(i+1)/np;
			f[1]++;
		}
		else if(pro[i]=='P'||pro[i]=='H')
		{
			pro1[i]=3;
			pos[2][f[2]]=(float)(i+1)/np;
			f[2]++;
		}
		else if(pro[i]=='A'||pro[i]=='G'||pro[i]=='Y'||pro[i]=='C'||pro[i]=='W')
		{
			pro1[i]=4;
			pos[3][f[3]]=(float)(i+1)/np;
			f[3]++;
		}
		else if(pro[i]=='L'||pro[i]=='V'||pro[i]=='I'||pro[i]=='F'||pro[i]=='M')
		{
			pro1[i]=5;
			pos[4][f[4]]=(float)(i+1)/np;
			f[4]++;
		}
		else
			printf("Error in data\n");
	}

//	fprintf(fout,">");
//	for(i=0;i<nn;i++)
//		fprintf(fout,"%c",name[i]);
//	stats(pro1,np,pro,name,nn);
	for(i=0;i<5;i++)
	{
		if(f[i]>max_freq)
			max_freq=f[i];
		for(j=0;j<f[i];j++)
			pos1[j]=pos[i][j];
		stats(pos1,f[i],m,order);
		fprintf(fout,"%lf ",(double)f[i]/MAX_freq);
		for(j=0;j<order;j++)
			fprintf(fout,"%e ",m[j]);
	}

	return;
}

void readfile(int argc, char *argv[])
{
	FILE *fin,*fout;
	char temp;
	char infile[50], outfile[50];
	char pro[M],name[N];	//array which stores the sequence of a protein
	int total=0;	//no of protein sequences
	int i=0,j=0;
	int ad=0,nonad=0;
	int pos[M];
	int nn=0,np=0;

	if(argc!=3)
	{
		if(argc!=2)
		{
			printf("Enter the inputfile: ");
			scanf("%s",infile);
		}
		else
			strcpy(infile,argv[1]);
		printf("Enter the outputfile: ");
		scanf("%s",outfile);
	}
	else
	{
		strcpy(infile,argv[1]);
		strcpy(outfile,argv[2]);
	}

	fin=fopen(infile,"r");
	if(fin==NULL)
	{
		printf("File %s not found.\n",infile);
		exit(3);
	}
	fout=fopen(outfile,"w");
	if(fout==NULL)
	{
		printf("Cannot open file %s to write.\n",outfile);
		exit(4);
	}

	total=initpos(fin,pos);

	for(j=0;j<total;j++)
	{
		fseek(fin,pos[j],SEEK_SET);
		for(i=0;fscanf(fin,"%c",&name[i])!=EOF&&name[i]!='\n';i++);
		nn=i;
		for(i=0;fscanf(fin,"%c",&pro[i])!=EOF&&pro[i]!='\n';i++);
		np=i;
		//fprintf(fout,"%lf ",(double)np/Lmax);
		procalc(fout,pro,np,name,nn);
		fprintf(fout,"%c\n",name[1]);
		if(name[1]=='1')
			ad++;
		else if(name[1]=='0')
			nonad++;
		else
		{
			printf("Error: In annotation 1 or 0 must be written before gi number.\n");
			exit(1);
		}
	}

	printf("Total no of proteins in the input = %d\n",total);
	printf("\t\t\tAdhesins\tNon-adhesins\n");
	printf("Total no\t:\t%d\t\t%d\n",ad,nonad);
	printf("Max freq = %d\n",max_freq);
	fclose(fin);
	fclose(fout);
	return;
}


int main(int argc, char *argv[])
{
	readfile(argc,argv);
	return 0;
}

