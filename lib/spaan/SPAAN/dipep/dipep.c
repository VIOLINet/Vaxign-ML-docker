#include<stdio.h>
#define M 10000
#define N 10000
#define cut_off 0.6
FILE *fin, *fout, *foutNN ;
char seq[M];
int STOP;
int adhesin;
float total_count[2][20][20];
int local_count[20][20];
int no_of_ad=0;
int no_of_nonad=0;
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



	//-------- adding the total dipeptide count ---------//
	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
		{
			total_count[adhesin][i][j] += (float)local_count[i][j]/(length-1);
		}


}











int main(int argc, char *argv[])
{
	int i,j,k,length;
	char temp, option, inputfile[50], outputfile[50], NNoutputfile[50];
	int significance[20][20];
	int significance_count=0;
	STOP=1;
	
	printf("Enter the input filename: ");
	scanf("%s",inputfile);
	printf("Enter the NN-outputfilename: ");
	scanf("%s",NNoutputfile);
	printf("Enter the filename for detailed output: ");
	scanf("%s",outputfile);
	fin=fopen(inputfile,"r");
	foutNN=fopen(NNoutputfile,"w");
	fout=fopen(outputfile,"w");

	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
		{
			significance[i][j]=0;
			for(k=0;k<2;k++)
				total_count[k][i][j]=0;
		}


	//----- counting the number of dipeps in the input file ------//
	fscanf(fin,"%c",&temp);
	while(STOP)
	{
		length=readfile(fin);
		dipep(length);
		if(adhesin==1)
			no_of_ad++;
		if(adhesin==0)
			no_of_nonad++;
	}



	//----- averaging over total no. of adhesins and nonadhesins ----//
	total_count[0][i][j] = total_count[0][i][j]/no_of_nonad;
	total_count[1][i][j] = total_count[1][i][j]/no_of_ad;



	//------- significance check --------//
	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
		{
			if((total_count[0][i][j] > total_count[1][i][j]))
			{
				if( ((total_count[0][i][j]-total_count[1][i][j])/total_count[0][i][j]) >= cut_off)
				{
					significance[i][j]=1;
					significance_count++;
				}
			}
			else
			{
				if((total_count[1][i][j] > total_count[0][i][j]))
					if( ((total_count[1][i][j]-total_count[0][i][j])/total_count[1][i][j]) >= cut_off)
					{
						significance[i][j]=1;
						significance_count++;
					}
			}
		}



	//---- debugging -----//
	printf("\n\n\n\n");
	for(i=0;i<20;i++)
		for(j=0;j<20;j++)
			if(significance[i][j]==1)
				printf("%c%c ",AA[i],AA[j]);

	printf("helu2");

	fseek(fin,0,SEEK_SET);
	STOP=1;
	fscanf(fin,"%c",&temp);
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




//	for(i=0;i<20;i++)
//		for(j=0;j<20;j++)
//			fprintf(foutNN,"count[%c][%c] = %f %f\n",AA[i],AA[j],(float)total_count[1][i][j]*100/no_of_ad,(float)total_count[0][i][j]*100/no_of_nonad);


	printf("NO. OF SIGNIFICANT DIPEPTIDES = %d\n",significance_count);
	fprintf(fout,"\n\n\n");
	printf("\n\n");
	fclose(fin);
	fclose(fout);
	fclose(foutNN);

	
	
	
	return 0;


}

