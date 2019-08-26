#include <stdio.h>
#define N 5

int main(int argc, char *argv[])
{
	FILE *fin[N], *fout, *fant;

	char antfile[50], infile1[50], infile2[50], infile3[50], infile4[50], infile5[50], infile6[50], outfile[50];
	int i, j, k;
	char c, c1, c2, s1[20], s2[20], s3[20];
	float p[6];
	int pfinal, ad, Tad=0, FP=0, TP=0;
	float cut_off, Sn=0, Sp=0, out;

//	fout=fopen(outfile,"w");
//	for(cut_off=0.00;cut_off<=1.0;cut_off+=0.01)
//	{
		TP=0;
		FP=0;
		Tad=0;
		for(k=0;k<N;k++)
			fin[k]=fopen(argv[k+1],"r");
		fant=fopen(argv[N+1],"r");
		fout=fopen(argv[N+2],"w");
		fprintf(fout,"SN\tPad-value\tProtein name (Annotation)\n");

//		printf("Enter the cut_off: ");
//		scanf("%f",&cut_off);

		i=0;
//		printf("%f\n",cut_off);
		while( fscanf(fin[0],"%s%s%s%f",s1,s2,s3,&p[0]) != EOF)
		{
			i++;
			for(k=1;k<N;k++)
				fscanf(fin[k],"%s%s%s%f",s1,s2,s3,&p[k]);
			out=p[0];
/*			for(j=1;j<6;j++)
				if(p[j]>=out)
					out=p[j];
*/			out=(p[0]*0.844444+p[1]*0.711111+p[2]*0.795556+p[3]*0.831111+p[4]*0.844444)/(0.844444+0.711111+0.795556+0.831111+0.844444);
/*			if(ad==1)
			{
				Tad++;
				if(out>=cut_off)
					TP++;
			}
			else
			{
				if(out>=cut_off)
					FP++;
			}
*/
			fprintf(fout,"%d\t",i);
//			for(k=0;k<N;k++)
//					fprintf(fout,"%f ",p[k]);
			fprintf(fout,"%f\t",out);
			do{
				fscanf(fant,"%c",&c);
				fprintf(fout,"%c",c);
			}while(c != '\n'); 
			do{
				fscanf(fant,"%c",&c);
			}while(c != '\n');
		}
//		Sn=(float)TP/Tad;
//		Sp=(float)TP/(TP+FP);
//		printf("FP=%d\n",FP);

//		fprintf(fout,"%f\t%f\t%f\n",cut_off,Sn,Sp);
		for(k=0;k<N;k++)
			fclose(fin[k]);
		fclose(fant);
//		fclose(fval);
//	}
	
	fclose(fout);

	return 0;
}
