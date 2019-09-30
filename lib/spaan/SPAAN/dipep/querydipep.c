#include<stdio.h>
#define M 15000
#define N 15000
FILE *fin, *fout, *foutNN, *finsign;
char seq[M];
int STOP;
int adhesin;
int local_count[20][20];
char AA[20]= {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};


int readfile(FILE *fin)
{
	char temp='>';
	int i=0,j=0,k1=1,k,length;


//	fprintf(fout,"\n\n");
//	printf("\n\n");
	if( temp == '>')
	{
		fscanf(fin,"%d",&adhesin);
		fprintf(fout,">");
		fprintf(fout,"%d",adhesin);
//		printf(">");
//		printf("%d",adhesin);
		do{
			fscanf(fin,"%c",&temp);
//			printf("%c",temp);
			fprintf(fout,"%c",temp);
		}while(temp!='\n');
	}



	while(k1)
	{
		if((fscanf(fin,"%c",&temp)) == EOF)
		{
			k1=0;
			STOP=0;
			break;
		}
		if(temp =='>')
		break;
		if(temp != '\n')
		{
			seq[i]=temp;
			i++;
		}
	 }
	return i;
}





void dipep(int length)
{
	char c1, c2;
	int i,j,k;
	int find1,find2;


	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
			local_count[i][j]=0;


	//-------- reading the file and counting the dipeptides --------//
	for(i=0;i<length-1;i++)
	{
		if(i==0)
			c1=seq[i];
		else
			c1=c2;
		c2=seq[i+1];



		//------- finding the dipeptide ---------//
		for(j=0;j<20;j++)
			if(AA[j]==c1)
				break;
		find1=j;
		for(j=0;j<20;j++)
			if(AA[j]==c2)
				break;
		find2=j;



		//--------- counting the dipeptides in the sequence ----------//
		local_count[find1][find2]++;
	}



}



int main(int argc,char *argv[])
{
	int i,j,k,length;
	char temp, option, inputfile[50], outputfile[50], NNoutputfile[50];
	int significance[20][20];
	STOP=1;
	
	if(argc!=3)
	{	printf("Enter the query filename: ");
		scanf("%s",inputfile);
		printf("Enter the NN-outputfilename: ");
		scanf("%s",NNoutputfile);
	}
	else
	{
		strcpy(inputfile,argv[1]);
		strcpy(NNoutputfile,argv[2]);
	}
	fin=fopen(inputfile,"r");
	foutNN=fopen(NNoutputfile,"w");
	fout=fopen("1.tmp","w");
	finsign=fopen("significant_dipep","r");

	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
		{
			significance[i][j]=0;
		}

	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
			fscanf(finsign,"%d",&significance[i][j]);



	while(STOP)
	{
		length=readfile(fin);
		dipep(length);
		for(i=0;i<20;i++)
			for(j=0;j<20;j++)
				if(significance[i][j] ==1)
					fprintf(foutNN,"%f ",(float)local_count[i][j]*100/(length-1));
		fprintf(foutNN," %d\n",adhesin);
	}

	fclose(fin);
	fclose(fout);
	fclose(foutNN);
	fclose(finsign);
	system("rm -f 1.tmp");
	return 0;

}

