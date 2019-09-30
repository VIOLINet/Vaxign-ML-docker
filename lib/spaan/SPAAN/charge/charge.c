#include<stdio.h>
#include<math.h>
#define M 15000
#define Lmax 15000	//expected maximum length of a protein
#define CLUmax 10	//expected maximum number of clusters in a protein
#define D 5
#define wmax np
#define wmin 50
#define gap 10

double fad=0,fnonad=0,fad_clu=0,fnonad_clu=0;
int ad=0,nonad=0;
int maxK=0;
//Pritilipi Adhikaar:
//			 Sachu, Bhootpurve chhatra,
//			 Bhartiya Prodyogiki Sansthan, Kanpur
//			 Dinaank: Prathma, Shukl Paksh,
//			 Asharh Mah, Vikram Samwat 2054
//			 Godhuli Vela, Saymkaal

void stats(float *s,int ns,double *m, int order)
{
	int i,j;
	m[0]=0;
	for(i=0;i<ns;i++)
		m[0]+=s[i];
	m[0]=m[0]/ns;
	for(j=1;j<order;j++)	//Its starting from 2nd order moment.
	{
		m[j]=0;
		for(i=0;i<ns;i++)
			m[j]+=pow(s[i]-m[0],j+1);
		m[j]=m[j]/ns;
	}
	return;
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
	fseek(fin,0,SEEK_SET);
	return total;
}

void procalc(FILE *fout,char *pro,int np,char *name,int nn,int type)
{
	char temp;
	int i=0,j=0;
	int pro1[M];
	float pos[M];
	int c=0,k=0,w,wprev;
	int clu_pos[100];
	int clu_len[100];
	int order=18;
	float f=0,f1=0,t,T,tprev;
	float clu_t[100];
	double f_clu=0;
	double m[20];


	for(i=0;i<np;i++)
	{
		pro1[i]=0;
		if(type==1)
		{
			if(pro[i]=='R'||pro[i]=='K'||pro[i]=='E'||pro[i]=='D')
			{
				pro1[i]=1;
				pos[c]=(float)(i+1)/np;
				c++;
			}
		}
		else if(type==2)
		{
			if(pro[i]=='R'||pro[i]=='K')
			{
				pro1[i]=1;
				pos[c]=(float)(i+1)/np;
				c++;
			}
		}
		else if(type==3)
		{
			if(pro[i]=='E'||pro[i]=='D')
			{
				pro1[i]=1;
				pos[c]=(float)(i+1)/np;
				c++;
			}
		}
		else
		{
			printf("Type mentioned is other than 1,2,3.\n");
			exit(1);
		}
	}
	f=(float)c/np;

	stats(pos,c,m,order);
	fprintf(fout,"%f ",f);
	for(i=0;i<order;i++)
		fprintf(fout,"%e ",m[i]);

	if(name[1]=='0')
	{
		fnonad+=f;
		fnonad_clu+=k;
		nonad++;
	}
	else if(name[1]=='1')
	{
		fad+=f;
		fad_clu+=k;
		ad++;
	}
	return;
}

void readfile(int argc, char *argv[])
{
	FILE *fin,*fout;
	char temp;
	char infile[50], outfile[50];
	char pro[M],name[200];	//array which stores the sequence of a protein
	int total=0;	//no of protein sequences
	int i=0,j=0;
	int type=1;
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

	//printf("Enter 1,2 or3 for selecting the type of charge clusters,\n\t1:\tmixed\n\t2:\tpositive\n\t3:\tnegative\n");
	//scanf("%d",&type);
	type=1;

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
		fprintf(fout,"%lf ",(double)np/Lmax);
		procalc(fout,pro,np,name,nn,type);
		fprintf(fout,"%c\n",name[1]);
	}

	printf("Total no of proteins in the input = %d\n",total);
	printf("\t\t\tAdhesins\tNon-adhesins\n");
	printf("Total no\t:\t%d\t\t%d\n",ad,nonad);
	printf("Avg Charge freq\t:\t%lf\t%lf\n",fad/ad,fnonad/nonad);
	fclose(fin);
	fclose(fout);
	return;
}

int main(int argc, char *argv[])
{
	readfile(argc,argv);
	return 0;
}
