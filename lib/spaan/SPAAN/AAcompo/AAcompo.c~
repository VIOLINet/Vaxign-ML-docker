#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define M 15000
#define D 5
#define wmax nn
#define wmin 1
#define gp1 10
#define SIZE 100
#define CLU 10
#define SHUFFLE 5000
#define STATS 4
#define MAX_freq 1000
FILE *fin,*fout,*fout1,*fout2;
int f_size[SIZE],n_clu[SIZE];
//Pritilipi Adhikaar:
//			 Sachu, Bhootpurve chhatra,
//			 Bhartiya Prodyogiki Sansthan, Kanpur
//			 Dinaank: Prathma, Shukl Paksh,
//			 Asharh Mah, Vikram Samwat 2054
//			 Godhuli Vela, Saymkaal

void freq(int *pro1,int np,char *name)
{
	int i;
	double f[20];
	for(i=0;i<20;i++)
		f[i]=0;
	for(i=0;i<np;i++)
		f[pro1[i]-1]++;
	for(i=0;i<20;i++)
		fprintf(fout,"%lf ",f[i]/MAX_freq);
	fprintf(fout,"%c\n",name[0]);
	return;
}

void stats(float *pro1,int np,char *name)
{
	int i,j;
	double mean=0,sd,mos,mok;
	double mm[STATS-1];
	for(i=0;i<np;i++)
		mean=mean+pro1[i];
	mean=mean/np;
//	fprintf(fout,"%lf ",mean);
	for(j=0;j<STATS-1;j++)
	{
		mm[j]=0;
		for(i=0;i<np;i++)
			mm[j]=mm[j]+pow(pro1[i]-mean,j+2);
		mm[j]=mm[j]/np;
//		fprintf(fout,"%lf ",mm[j]);
	}
	sd=sqrt(mm[0]);
	mos=mm[1]/pow(mm[0],1.5);
	mok=mm[2]/pow(mm[0],2);
		 
	fprintf(fout,"%lf %lf %lf %lf %c\n",mean,mm[0],mm[1],mm[2],name[0]);
	return;
}

void convert(char *pro,int np,char *name,int nn)
{
	int i,j;
	int pro1[M];
	for(i=0;i<np;i++)
	{
		if(pro[i]=='Q')
			pro1[i]=1;
		else if(pro[i]=='E')
			pro1[i]=2;
		else if(pro[i]=='D')
			pro1[i]=3;
		else if(pro[i]=='R')
			pro1[i]=4;
		else if(pro[i]=='K')
			pro1[i]=5;
		else if(pro[i]=='H')
     			pro1[i]=6;
		else if(pro[i]=='Y')
			pro1[i]=7;
		else if(pro[i]=='W')
			pro1[i]=8;
		else if(pro[i]=='F')
			pro1[i]=9;
		else if(pro[i]=='M')
			pro1[i]=10;
		else if(pro[i]=='L')
			pro1[i]=11;
		else if(pro[i]=='I')
			pro1[i]=12;
		else if(pro[i]=='V')
			pro1[i]=13;
		else if(pro[i]=='C')
			pro1[i]=14;
		else if(pro[i]=='P')
			pro1[i]=15;
		else if(pro[i]=='A')
			pro1[i]=16;
		else if(pro[i]=='G')
			pro1[i]=17;
		else if(pro[i]=='T')
			pro1[i]=18;
		else if(pro[i]=='S')
			pro1[i]=19;
		else if(pro[i]=='N')
			pro1[i]=20;
		else
		{
			printf("Error: There is an undefined amino acid in the given protein sequence\n");
			for(j=0;j<nn;j++)
				printf("%c",name[j]);
			return;
		}
//		pro1[i]=pro1[i]/10;
	}
	freq(pro1,np,name);
	return;
}

void shuffle(char *pro,int np)
{
	int i,j,r;
	char temp;
	for(i=0;i<SHUFFLE;i++)
	{
		srand(time(NULL));
		for(j=0;j<np;j++)
		{
			r =(float)(np-1-j)*rand()/RAND_MAX;
//			printf("%d ",r);
			r=r+j;
			temp=pro[j];
			pro[j]=pro[r];
			pro[r]=temp;
		}	
	}
	return ;
}

void freq_size(int nclu,int *clu_len)
{
	int i,j;
	for(i=0;i<nclu;i++)
		if(clu_len[i]<=SIZE)
			f_size[clu_len[i]-1]++;
	if(nclu<=CLU)
		n_clu[nclu]++;
	return;
}
void procalc(char *pro,int np,char *name,int nn,int type)
{	
	int i=0,j=0;	
	int pro1[M];
	int c=0,k=0,w,wprev;
	float f=0,f1=0,t,T,tprev;
	double f_clu=0;
	int clu_pos[100];
	int clu_len[100];
	float clu_t[100];
		
/*	for(i=0;i<np;i++)
		printf("%c",pro[i]);
	printf("\n");
*/	
	shuffle(pro,np);
	
/*	for(i=0;i<np;i++)
		printf("%c",pro[i]);
	printf("\n");
*/
	for(i=0;i<np;i++)
	{
		pro1[i]=0;
		if(type==1)
		{
			if(pro[i]=='R'||pro[i]=='K'||pro[i]=='E'||pro[i]=='D')
			{
				pro1[i]=1;
				c++;
			}
		}
		else if(type==2)
		{
			if(pro[i]=='R'||pro[i]=='K')
			{
				pro1[i]=1;
				c++;
			}
		}
		else if(type==3)
		{
			if(pro[i]=='E'||pro[i]=='D')
			{
				pro1[i]=1;
				c++;
			}
		}
	}
	f=(float)c/np;
	
	if(np>=1500) T=5;
	else if(np<1500 && np>=750) T=4.5;
	else if(np<750) T=4;

	for(i=0;i<np;i++)
	{
		if(pro1[i]==1)
		{
			c=1;
			wprev=1;
			tprev=0;
			clu_len[k]=0;
			for(w=2;w<=wmax && i+w<np;w++)
			{
				if(pro1[i+w-1]==1)
				{
				//	if(w-wprev>gp1)
				//		break;
					c++;
					wprev=w;						
					if(w>=wmin)
					{
						t=(c-w*f)/sqrt(w*f*(1-f));
						if(t>T)
							if(t>=tprev)
							{
								clu_len[k]=w;
								clu_t[k]=t;
								tprev=t;
							}
					}
				}
			}
			if(clu_len[k]>=wmin)
			{
				if(clu_t[k]>=T)
				{
					if(k==0)
						clu_pos[k++]=i+1;
					else
					if(i+1>clu_pos[k-1]+clu_len[k-1]-1)
						clu_pos[k++]=i+1;
					else
						if(clu_t[k]>clu_t[k-1])
						{
							clu_pos[k-1]=i+1;
							clu_len[k-1]=clu_len[k];
							clu_t[k-1]=clu_t[k];
						}
				}
			}
		}
	}

	freq_size(k,clu_len);	
/*	for(i=0;i<k;i++)
		f_clu = f_clu + clu_len[i];
	f_clu=f_clu/np;
	
	fprintf(fout,"%f\t%lf\n",f,f_clu);
*/
	return;
}

void readfile(int argc, char *argv[])
{
	char temp;
	char infile[20];
	char outfile[20];
	char pro[M],name[200];	//array which stores the sequence of a protein
	int n=0;	//no of protein sequences
	int i=0,j=0;
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
//	fout1=fopen("out1.txt","w");

	for(i=0;i<SIZE;i++)
		f_size[i]=0;
	for(i=0;i<CLU;i++)
		n_clu[i]=0;
	i=0;
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
				convert(pro,i,name,j);
			i=0;
			n++;
			j=0;
			do{
				fscanf(fin,"%c",&temp);
				name[j++]=temp;
			}while(temp!='\n');
		}
	}
	convert(pro,i,name,j);
	printf("Total no of proteins in the input: %d\n",n);
/*	for(i=0;i<SIZE;i++)
		fprintf(fout,"%d %d\n",i+1,f_size[i]);
	for(i=0;i<CLU;i++)
		fprintf(fout1,"%d %d\n",i,n_clu[i]);
*/
	fclose(fin);
	fclose(fout);
//	fclose(fout1);
	return;
}

int main(int argc, char *argv[])
{
	readfile(argc,argv);
	return 0;
}
