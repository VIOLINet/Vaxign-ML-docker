#include<stdio.h>
#include<math.h>
// #include<stdlib.h>
#define M 15000
#define N 15000
#define iteration_limit 1000
// #define count_limit 5
FILE *fin, *fout, *NNinput;
char seq[M];
int STOP;
int rand_seed=2;
int count_limit;
int length_limit;
float total_freq[2][20];
float total_frac[2][20];
int number_of_adhesins = 0;
int number_of_nonadhesins =0;
char AA[20]= {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
int adhesin=1;

int readfile(FILE *fin)
{
	char temp='>';
	int i=0,j=0,k1=1,k,length;


	//fprintf(fout,"\n\n");
//	printf("\n\n");
	if( temp == '>')
	{
		fscanf(fin,"%d",&adhesin);
//		fprintf(fout,">");
		//fprintf(fout,"\n\n");
		do{
			fscanf(fin,"%c",&temp);
//			printf("%c",temp);
			//fprintf(fout,"%c",temp);
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





int rand1(int length)
// from K&R
//- produces a random number between 0 and length.
{
rand_seed = rand_seed * 1103515245 +12345;
return (unsigned int)(rand_seed / 65536) % (length+1);
}




void repeat(int length,FILE *NNinput)
{
	int i, j, k, tmp, option;
	int start,stop;
	int count_of_repeats=0, repeat_length, decision;
	float fraction=0.0;
	int doublet=1, length_of_repeat[20], repeat_count[20];
	float freq[2][20];
	float frac[2][20];


	for(i=0;i<20;i++)
	{
		length_of_repeat[i]=0;
		repeat_count[i]=0;
		freq[0][i]=0.0;	freq[1][i]=0.0;
		frac[0][i]=0.0;	frac[1][i]=0.0;
	}


		//--------- finding the repeats ------------//
	for(i=0;i<length;)
	{
		if(count_of_repeats==0)
		{
			tmp=i;
			decision=1;
		}
		if(count_of_repeats==1)
		{
//			if(tmp<(i+length_limit))
//				tmp++;
			if(tmp>=(i+length_limit))
			{
				i++;
				tmp=i;
				decision=1;
			}
		}
	
		
		do{
				tmp++;
				if(tmp>=length)
				{
					decision=0;
					break;
				}
		   }while(seq[tmp] != seq[i]);
		repeat_length=tmp-i;           // if repeat is ACDACDACD... then repeat_length=3 
		count_of_repeats=1;


			//------------ counting the number of repeats -------------//
		while(decision)
		{
			if(repeat_length==1)
			{
				do{
						count_of_repeats++;
				}while(seq[i]==seq[i+count_of_repeats]);
				decision=0;
			}
			else
			{
				for(j=0;j<repeat_length;j++)
					if(seq[i+j] != seq[i+j+(count_of_repeats*repeat_length)])
						decision=0;
			}
			if(decision==1)
				count_of_repeats++;
 		}

		
		
		
		if(count_of_repeats >= count_limit)
		{
			start =i;
			stop = i+ (count_of_repeats*repeat_length) -1;
		}

		//----------- printing the results -------------//
		if(count_of_repeats >= count_limit)
		{
			fraction += (float)(stop-start+1);
//			printf("repeat found at position %d of length %d & count= %d ",i,repeat_length,count_of_repeats);
			//fprintf(fout,"repeat found at position %d of length %d & count= %d ",i,repeat_length,count_of_repeats);
			for(j=start;j<=stop;j++)
			{
//				printf("%c",seq[j]);
				//fprintf(fout,"%c",seq[j]);
			}
//			printf("\n");
//			fprintf(fout,"\n");
		}
	

	
		//--------- calculating the frequency ---------//
		if(count_of_repeats >= count_limit)
		{
			doublet=1;
			for(j=start;j<=stop;j++)
				if(seq[j] != seq[start])
					doublet=0;
			if(doublet==1)
			{
				for(j=0;j<20;j++)
				{
					if(AA[j]==seq[start])
						break;
				}
				length_of_repeat[j] += (stop-start+1);
				repeat_count[j]++;
				//fprintf(fout,"\nlength_of_repeat[%c]=%d\n",AA[j],length_of_repeat[j]);
			}
		}

	
		if(count_of_repeats>1)
		{
			i += ((count_of_repeats-1)*repeat_length) ;	
			tmp=i;
			decision=1;
		}
	}
	fraction = fraction/length;
	//fprintf(fout,"fraction = %.2f\n",fraction);
	for(j=0;j<20;j++)
	{
		frac[adhesin][j] += ((float)length_of_repeat[j]/length);
		freq[adhesin][j] += ((float)repeat_count[j]/length);
		total_frac[adhesin][j] += ((float)length_of_repeat[j]/length);
		total_freq[adhesin][j] += ((float)repeat_count[j]/length);

//		fprintf(fout,"length_of_repeat[%c]=%d repeat_count[%c]=%d frac[%c]=%.2f ",AA[j],length_of_repeat[j],AA[j],repeat_count[j],AA[j],frac[j]);
		//fprintf(fout," frac[%d][%c]=%f ",adhesin,AA[j],frac[adhesin][j]);
		fprintf(NNinput,"%f ",100.0*frac[adhesin][j]);
	}
	fprintf(NNinput,"%d\n",adhesin);
	//fprintf(fout,"\n");





}













int main(int argc, char *argv[])
{
	int i,j,length,k,choice;
	char temp,option, infile[50], outfile[50], NNinputfile[50];
//	fin=fopen("twogene.txt","r");
//	fout=fopen("out.txt","w");
	STOP=1;
//	printf("Hi1\n");
	for(i=0;i<2;i++)
		for(j=0;j<20;j++)
		{
			total_freq[i][j]=0.0;
			total_frac[i][j]=0.0;
		}
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
	//printf("Enter the outputfilename: ");
	//scanf("%s",outputfile);

	fin=fopen(infile,"r");
	if(fin==NULL)
	{
		printf("File %s not found.\n",infile);
		exit(3);
	}

	NNinput=fopen(outfile,"w");
	if(NNinput==NULL)
	{
		printf("Cannot open file %s to write.\n",outfile);
		exit(4);
	}
	//fout=fopen(outputfile,"w");
	count_limit=2;
	length_limit=10;
	//system("clear");

start:

//	printf("\t\tDEFAULT count_limit = 2 & length_limit = 10\n");
//	printf("\t\t 1. Accept 2. change the default values\n");
//	printf("\t\t Enter the option :");
//	scanf("%d",&choice);
	choice=1;
	if(choice==2)
	{
		printf("Enter the lower limit for repeat_count: ");
		scanf("%d",&count_limit);
		printf("Enter the lower limit for repeat_length: ");
		scanf("%d",&length_limit);
	}
//	if(choice==1)
//		printf("\t\tdefault values ACCEPTED\n");
	if(choice != 1 && choice != 2)
	{
		printf("\t\tincorrect option\n");
		goto start;
	}
	fscanf(fin,"%c",&temp);
	while(STOP)
	{
		length=readfile(fin);
		//fprintf(fout,"length of the protein: %d\n",length);
 		repeat(length,NNinput);
		if(adhesin==0)
			number_of_nonadhesins++;
		if(adhesin==1)
			number_of_adhesins++;
//		printf("%d ",adhesin);
//		fprintf(fout,"\n\n\tWISH TO GO TO NEXT PROTEIN: (Y/N) ");
//		scanf("%c",&option);
//		if(option=='n')
//			STOP=0;
//
//		else
//			continue;
	}
//	printf("\n");
	for(j=0;j<20;j++)
	{
		total_frac[0][j] = ((float)total_frac[0][j]/number_of_nonadhesins);
		total_frac[1][j] = ((float)total_frac[1][j]/number_of_adhesins);
		total_freq[0][j] = ((float)total_freq[0][j]/number_of_nonadhesins);
		total_freq[1][j] = ((float)total_freq[1][j]/number_of_adhesins);
//		printf("frac[0][%c]=%f frac[1][%c]=%f	freq[0][%c]=%f freq[1][%c]=%f\n",AA[j],total_frac[0][j]*100.0,AA[j],total_frac[1][j]*100.0,AA[j],total_freq[0][j]*100.0,AA[j],total_freq[1][j]*100.0);
		//fprintf(fout,"frac[0][%c]=%f frac[1][%c]=%f	freq[0][%c]=%f freq[1][%c]=%f\n",AA[j],total_frac[0][j]*100.0,AA[j],total_frac[1][j]*100.0,AA[j],total_freq[0][j]*100.0,AA[j],total_freq[1][j]*100.0);

	}

	//fprintf(fout,"\n\n\n");
//	printf("\n\n");
	fclose(fin);
    	//fclose(fout);
	fclose(NNinput);
	return 0;
}































