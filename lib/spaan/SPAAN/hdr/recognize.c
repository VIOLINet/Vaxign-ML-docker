#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>		/* for setting random number seed */
#include <unistd.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define min3(a,b,c) (min((min((a), (b))), (c)))

#define FREE(ptr) if (ptr!=NULL) free(ptr);
#define MAXFILES 1000

#define NUMCATEGORIES 8 // Number of different phoneme categories (# of outputs of NN)
#define MAXWINDOWS 300000 // Max number of windows in a file to recognize
#define MAXTAGS 10 // Max length of a word's phonetic description
#define INF 100000 // Infinity !!

#define DEBUG 0 // Used in parts of the program written by RR

#define PH_SILENCE 0

int debug= 1;

struct pattern_struct {
  float *data;
  struct pattern_struct *next;
};
typedef struct pattern_struct Pattern;
Pattern *firstPattern;
Pattern *lastPattern;

struct set_struct {
  int num;
  Pattern **patterns;
/*  float *targets;*/
};
typedef struct set_struct Set;

int linearOut=0;		/* 1 if output is linear. */
int InputsThrough=0;		/* 1 to pass inputs through to last layer */

int numInputs, numOutputs, numHidden1, numHidden2=0, numOutputInputs;
int epochs, epochPrintInterval=0;
float correctCrit;

float *wh1;		/* firts hidden layer units' weights */
float *wh2;		/* optional second layer hidden units' weights */
float *wo;		/* output unit's weights */
float *dwh1;
float *dwh2;
float *dwo;
float *wp, *swp;
float *x;			/* input vector */
float *yh1;			/* hidden units' outputs */
float *yh2;			/* second layer hidden units' outputs */
float *yo;			/* output unit's output  */
//float *deltao, *deltah1, *deltah2, *deltax; /* derivatives */
//float *whchange = NULL;
//float *wochange = NULL;
//float *error;
//float epoch_error;
//float oRate, hRate, mom;
//float *swh1, *swh2, *swo;	/* best weights so far */
//float lowestValidateError;	/* keep track of best so far */
//int bestEpoch;

int highPrecision= 0;		/* 1 to use 32 bit backprop */

/* Remember all files and associated binary versions so later experiments */
/* can use them without converting again. */
struct file_struct {
  char name[100];
  char binName[100];		/* unused */
  int num;
  char binaryName[100];
  Pattern *first;
};
typedef struct file_struct DataFile;

DataFile dataFiles[MAXFILES];
int numDataFiles;
FILE *specFile;

int phone_to_plot;
char *melFile, *netFile, *outputFile, *graphFile; //, *validateFiles[100], *testFiles[100];
int nTrain, nValidate=0, nTest;
/*char resultsFile[100] = " ";*/
int summarize=0;		/* if 1, print one line per run */

Set mel; //, testing, validating;

char previousAlgorithm[100] = "none";
char newAlgorithm[100];
char *commandName;

int numwindows;
float prob[MAXWINDOWS][NUMCATEGORIES];

typedef struct {
  int numtags;
  int tag[MAXTAGS]; // Max no of phones in a word is MAXTAGS
  char word[10];
} Template;

Template T[20];
int numtemplates;

/*** function prototypes ***/

int fexists(char *filename);
int mysystem(char *s);		/* delay after doing command */
void initParse(char *file);
int parseNext(void);
int findFile(char *filename);
void addIfNew(char *filename);
void usage(void);
void loadSet(Set *set, char *file);
void readPatterns(DataFile *dataFile);

void dumpOutput(Set *mel, char *outputFile);

void calcOutput(float *x, float *wh1, float *wh2, float *wo);
//void printNet(int ni, int nh1, int nh2, int no, float *x, float *wh1,
//      float *wh2, float *wo, float *yh1, float *yh2, float *yo);
float logistic(float x);
float Dlogistic(float y);
void nLinesWords (char *filename, int *rows, int *words);
void init_acgrnd_by_time();
float acgrnd();
char *date(void);

/** This code uses NN evaluation functions from code by
    Charles W. Anderson, 1996
    Provides specification file structure
**/

/**********************************************************************/

float *float_malloc(int nb) {
  void *p;
  float *q;
  p = malloc(nb*4);
  q = (float *)p;
  return q;
}

void bye(char *msg) {
  printf("ERROR: %s\n", msg);
  exit(0);
}

/**********************************************************************/

void init_templates() {
  numtemplates = 3;

  strcpy(T[0].word, "ek");
  T[0].tag[0]=0;
  T[0].tag[1]=1;
  T[0].tag[2]=0;
  T[0].tag[3]=2;
  T[0].tag[4]=0;
  T[0].numtags=5;

  strcpy(T[1].word, "do");
  T[1].tag[0]=0;
  T[1].tag[1]=3;
  T[1].tag[2]=4;
  T[1].tag[3]=0;
  T[1].numtags=4;

  strcpy(T[2].word, "teen");
  T[2].tag[0]=0;
  T[2].tag[1]=5;
  T[2].tag[2]=6;
  T[2].tag[3]=7;
  T[2].tag[4]=0;
  T[2].numtags=5;
}

/**********************************************************************/
// Use current prob[][] array in both these

// Local phone-phone distance between computed probabilites prob[i] and actual T[tnum].tag[j]
float d(int tnum, int i, int j) {
  int k, rank;

  // Find rank of T[tnum][j] in prob[i] (0..NUMCATEGORIES-1)
  rank=1;
  for (k=0; k<NUMCATEGORIES; k++) {
    if (k==T[tnum].tag[j]) continue;
    if (prob[i][k]>prob[i][T[tnum].tag[j]])
      rank++;
  }

  // Penalise for too much silence (should not be needed with balanced training sets)
  if (T[tnum].tag[j]==PH_SILENCE) {
    if (rank>1) return 2.25*rank;
  }
  return rank;
}

// Global distance from template... Uses Dynamic Time Warping //
float distance_from_template(int tnum) {
  float D[MAXWINDOWS][MAXTAGS];
  int p[MAXWINDOWS][MAXTAGS];
  int q[MAXWINDOWS][MAXTAGS];
  int i, j, safei, prevq;
  float dx, dy;
  float p01, p10, p11; 

  dx=0.2;
  dy=0.2;

  D[0][0]=0;
  p[0][0]=-1;
  q[0][0]=-1;
  for (i=1; i<T[tnum].numtags; i++) {
    D[0][i]=INF;
    p[0][i]=-1;
    q[0][i]=-1;
  }

  for (i=1; i<numwindows; i++) {
    D[i][0] = D[i-1][0] + d(tnum, i, 0);
    p[i][0]=i-1;
    q[i][0]=0;
    for (j=1; j<T[tnum].numtags; j++) {
      p11 = D[i-1][j-1] + 2 * d(tnum, i, j);
      p10 = D[i-1][j] + d(tnum, i, j) + dx;
      p01 = D[i][j-1] + d(tnum, i, j) + dy;
      D[i][j] = min3(p11, p10, p01);
      if (D[i][j]==p11) {
	p[i][j]=i-1;
	q[i][j]=j-1;
      }
      else if (D[i][j]==p10) {
	p[i][j]=i-1;
	q[i][j]=j;
      }
      else if (D[i][j]==p01) {
	p[i][j]=i;
	q[i][j]=j-1;
      }
      else {
	bye("Bad bug in distance_from_template");
      }
    }
  }

  printf("Template %s, %d\n", T[tnum].word, tnum);
  i=numwindows-1; 
  j=T[tnum].numtags-1;
  prevq = q[i][j];
  printf("(%d, %d)  ", i, j);
  while (i!=-1 || j!=-1) {
    if ((q[i][j]!=prevq) || (DEBUG)) {
      printf("(%d, %d)  ", i, j);
      prevq = j;
    }
    safei=p[i][j];
    j=q[i][j];
    i=safei;
  }
  printf("\n");

  return D[numwindows-1][T[tnum].numtags-1];
}

/**********************************************************************/
int main(int argc, char *argv[]) {
  int u, i; // , j, p, pattern;
  //int ep;
  //int junk;
  //float sum, junkFloat;
  FILE *netfile, *graphfile;
  int oldni, oldnh1, oldnh2, oldno;
  char string[200];
  float min_dist, min_dist2, curr_dist;
  int min_dist_template;


  commandName = argv[0];

  if (argc != 2) {
    usage();
    exit(1);
  }

  /* Log file */
  /*
    if ((filep=fopen("train.log","a")) != NULL) {
    for (i=0; i<argc; i++) 
      fprintf(filep,"%s ", argv[i]);
    fprintf(filep,"\n");
    fclose(filep);
    }
  */

  //init_acgrnd_by_time();

  fexists(argv[1]);
  initParse(argv[1]);

  while (parseNext()) {

    numOutputInputs = max(numHidden1,numHidden2);

    if (numOutputInputs == 0 || InputsThrough)
      numOutputInputs += numInputs;

    /* Load patterns */

    loadSet(&mel,melFile);

    /* Set up other data structures */

    if (oldni != numInputs || oldnh1 != numHidden1 || 
	oldnh2 != numHidden2 || oldno != numOutputs) {

      FREE(yo); FREE(yh1);  FREE(wh1); FREE(wo);
      FREE(dwh1); FREE(dwo);
      if (yh2 != NULL) {
      	FREE(yh2); FREE(wh2); FREE(dwh2);
      }


      /*
      yo = (float *) malloc (numOutputs * sizeof(float));
      yh1 = (float *) malloc (numHidden1 * sizeof(float));
      error = (float *) malloc (numOutputs * sizeof(float));
      wh1 = (float *) malloc((numInputs+1)*numHidden1 * sizeof(float));
      wo = (float *) malloc(numOutputs*(numOutputInputs+1) * sizeof(float));
      dwh1 = (float *) malloc((numInputs+1)*numHidden1 * sizeof(float));
      dwo = (float *) malloc(numOutputs*(numOutputInputs+1) * sizeof(float));
      swh1 = (float *) malloc((numInputs+1)*numHidden1 * sizeof(float));
      swo = (float *) malloc(numOutputs*(numOutputInputs+1) * sizeof(float));
      if (numHidden2 > 0) {
	yh2 = (float *) malloc (numHidden2 * sizeof(float));
	wh2 = (float *) malloc((numHidden1+1)*numHidden2 * sizeof(float));
	dwh2 = (float *) malloc((numHidden1+1)*numHidden2 * sizeof(float));
	swh2 = (float *) malloc((numHidden1+1)*numHidden2 * sizeof(float));
      } else {
	yh2 = wh2 = dwh2 = swh2 = NULL;
      }
      */

      yo = float_malloc (numOutputs * sizeof(float));
      yh1 = float_malloc (numHidden1 * sizeof(float));
      wh1 = float_malloc((numInputs+1)*numHidden1 * sizeof(float));
      wo = float_malloc(numOutputs*(numOutputInputs+1) * sizeof(float));
      dwh1 = float_malloc((numInputs+1)*numHidden1 * sizeof(float));
      dwo = float_malloc(numOutputs*(numOutputInputs+1) * sizeof(float));
      if (numHidden2 > 0) {
	yh2 = float_malloc (numHidden2 * sizeof(float));
	wh2 = float_malloc((numHidden1+1)*numHidden2 * sizeof(float));
	dwh2 = float_malloc((numHidden1+1)*numHidden2 * sizeof(float));
      } else {
	yh2 = wh2 = dwh2 = NULL;
      }

    }
    oldni = numInputs;
    oldnh1 = numHidden1;
    oldnh2 = numHidden2;
    oldno = numOutputs;

    /* Initialize previous weight changes to zero. */


    //printf("\nWeights\n");
    if ((netfile = fopen(netFile, "r"))==NULL) {
      printf("Couldn't open %s\n", netFile);
      exit(1);
    }
    while (fscanf(netfile, "%s", string) && strncmp(string, "Weights", 7)!=0);
    
    for (u=0; u<numHidden1; u++) {
      for (i=0; i<numInputs+1; i++)
	fscanf(netfile, "%f",(wh1+u*(numInputs+1)+i));
      //printf("%f\n",*(swh1+u*(numInputs+1)+i));
      //printf("\n");
    }

    if (numHidden2 > 0) {
      //printf("\n");
      for (u=0; u<numHidden2; u++) {
	for (i=0; i<numHidden1+1; i++)
	  fscanf(netfile, "%f",(wh2+u*(numHidden1+1)+i));
	//printf("\n");
      }
    }

    //printf("\n");
    for (u=0; u<numOutputs; u++) {
      for (i=0; i<numOutputInputs+1; i++)
	fscanf(netfile, "%f",(wo+u*(numOutputInputs+1)+i));
      //printf("\n");
    }

    break;
  }

  //init_templates();

  dumpOutput(&mel, outputFile);

  /***********
  // To draw graphs of a phone's probability against time
  if (strcmp(graphFile, "all")!=0) {
    printf("Writing phone probability plot file %s...\n", graphFile);
    graphfile = fopen(graphFile, "w");
    if (graphfile==NULL) bye("Couldn't open graph file");
    for (i=0; i<min(50, numwindows); i++) {
      fprintf(graphfile, "%d %.3f\n", i, prob[i][phone_to_plot]);
    }
    fclose(graphfile);
  }
  else {
    printf("Writing ALL phone probability plot files...\n");
    for (u=0; u<NUMCATEGORIES; u++) { 
      sprintf(string, "%d.prob", u);
      graphfile = fopen(string, "w");
      if (graphfile==NULL) bye("Couldn't open graph file");
      for (i=0; i<min(100, numwindows); i++) {
	fprintf(graphfile, "%d %.3f\n", i, prob[i][u]);
      }
      fclose(graphfile);
    }
  }

  //numwindows=80;
  // Find closest template by Dynamic Time Warping
  min_dist=INF;
  min_dist2=INF;
  for (i=0; i<numtemplates; i++) {
    curr_dist = distance_from_template(i);
    if (curr_dist<min_dist) {
      min_dist2=min_dist;
      min_dist=curr_dist;
      min_dist_template=i;
    }
    else if (curr_dist<min_dist2) {
      min_dist2=curr_dist;
    }
    printf("%s : %.2f\n", T[i].word, curr_dist);
  }


  // Estimate the confidence in the result,
  // and then decide if the output is reliable enough to show to user
  if ((min_dist<1.75*numwindows && min_dist2-min_dist>numwindows/4)|| 
      ((min_dist2-min_dist>numwindows/2) && (min_dist<2.5*numwindows))) {
    printf("\nOUTPUT : %s\n\n", T[min_dist_template].word);
  }
  else {
    printf("\nOUTPUT : Unable to classify\n\n");
  }
  *********/

  return 0;
}


/**********************************************************************/

#define IF(str) if (strncmp(token,str,strlen(str))==0)
#define EIF(str) else if (strncmp(token,str,strlen(str))==0)
#define NEXTTOKEN(what) if (fscanf(specFile,"%s",token) == EOF) { printf("Specification file ended while %s\n",what); exit(1);}
#define NEXTF(next) atof(token)
#define NEXTI(next) atoi(token)


void initParse(char *file) {
  if ((specFile=fopen(file,"r")) == NULL) {
    usage();
    exit(1);
  }
  numDataFiles = 0;
  phone_to_plot=0;
}

/**********************************************************************/

int parseNext(void) {
  static char token[1000];
  int notEOF = 1, advance;

  notEOF = fscanf(specFile,"%s",token) != EOF;

  while(notEOF) {

/*    printf("Token is %s\n",token);*/

    advance = 1;
    IF("-dat") {
      NEXTTOKEN("reading dat file name");
      melFile = (char *) malloc(sizeof(char)*(1+strlen(token)));
      strcpy(melFile,token);
      addIfNew(melFile);
    }
    EIF("-net") { 
      NEXTTOKEN("reading net file name");
      netFile = (char *) malloc(sizeof(char)*(1+strlen(token)));
      strcpy(netFile,token);
    }
    EIF("-output_file") { 
      NEXTTOKEN("reading output file name");
      outputFile = (char *) malloc(sizeof(char)*(1+strlen(token)));
      strcpy(outputFile,token);
    }
    EIF("-graph_file") { 
      NEXTTOKEN("phone to plot");
      phone_to_plot = atoi(token);
      NEXTTOKEN("reading graph file name");
      graphFile = (char *) malloc(sizeof(char)*(1+strlen(token)));
      strcpy(graphFile,token);
    }
    EIF("-nhiddens1") {NEXTTOKEN("number of hidden units 1"); 
		       numHidden1 = atoi(token);}
    EIF("-nhiddens2") {NEXTTOKEN("number of hidden units 2"); 
		       numHidden2 = atoi(token);}
    EIF("-ni") {NEXTTOKEN("number of inputs"); numInputs = atoi(token);}
    EIF("-no") {
      NEXTTOKEN("number of outputs"); numOutputs = atoi(token);
      //if (numOutputs!=NUMCATEGORIES) bye("Incompatible # of outputs and phoneme types");
    }
/**    EIF("-re") {NEXTTOKEN("result file name"); strcpy(resultsFile,token);}**/
    EIF("-end") {return 1;}	/* signal that experiment may start */
    else if (token[0] == '-') { /* argument has - but didn't match */
      printf("Unrecognized argument %s\n",token);
      usage();
      exit(1);
    }
    if (advance) 
      notEOF = fscanf(specFile,"%s",token) != EOF;
  }
  return 0;			/* signal that no more tokens */
}
/**********************************************************************/
int findFile(char *filename) {
  int i; //, found=0;
  
  for (i=0; i<numDataFiles; i++)
    if (strcmp(dataFiles[i].name,filename) == 0) {
      return i;
    }

  return -1;
}

/**********************************************************************
 * Match filename against already read files.  If not found, add it.
 **********************************************************************/
void addIfNew(char *filename) {
  int i, lines, words;
  //char string[200];

  if (debug) {printf("Reading %s\n",filename); fflush(stdout);}

  if((i = findFile(filename)) == -1) {
    if (numDataFiles >= MAXFILES) {
      printf("Already have the maximum of %d data files.\n",MAXFILES);
      exit(1);
    }
    strcpy(dataFiles[numDataFiles].name,filename);
    nLinesWords(filename,&lines,&words);
    dataFiles[numDataFiles].num = words / (float)(numInputs+numOutputs);
    if (dataFiles[numDataFiles].num*(numInputs+(float)numOutputs) != (float)words) {
      printf("For file %s, number of words not divisible by num inputs plus \
num outputs.\n",dataFiles[numDataFiles].name);
      printf(" %d Words  %d Inputs %d Outputs %d rows\n",
	     words,numInputs,numOutputs,dataFiles[numDataFiles].num);
      exit(1);
    }
    if (debug) {
      printf("%s l %d w %d num %d\n",dataFiles[numDataFiles].name,lines,words,
       dataFiles[numDataFiles].num);
      fflush(stdout);
    }

    readPatterns(&dataFiles[numDataFiles]);

    numDataFiles++;
  }
}

/**********************************************************************/

void usage(void) {
  printf("Usage: recognize <spec-file>\n");
  //<mel-to-recognize> <NN-weights-file>\n");
}

/**********************************************************************
 * Place pointers to all patterns into the training, validation, or
 * testing set.
 **********************************************************************/
void loadSet(Set *set, char *file) {
  int i, j, np,  nFiles;
  //int ret, jj;
  Pattern *p;
  int *fileIndices;
  nFiles=1;
  fileIndices = (int *) malloc(nFiles*sizeof(int));
  //char dest[200], string[200];

  set->num = 0;
  for (i=0; i<nFiles; i++) {
    fileIndices[i] = findFile(file);
    /*printf("findfile %s index %d\n",files[i],fileIndices[i]);*/
    set->num += dataFiles[fileIndices[i]].num;
  }

  //if (set->num>0)
  set->patterns = (Pattern **) malloc(set->num * sizeof(Pattern *));

  np = 0;
  for (i=0; i<nFiles; i++) {
    p = dataFiles[fileIndices[i]].first;
    for (j=0; j<dataFiles[fileIndices[i]].num; j++){
      set->patterns[np++] = p;
      p = p->next;
    }
  }

/**
  printf("Printing set\n");
  for (i=0; i<set->num; i++) {
    printf("Pattern %d ",i);
    printData(" ",set->patterns[i]);
    printf("\n");
  }
**/

}

int printData(char *s, Pattern *p) { /*"added pattern",first+j);*/
  int i;
  printf("%s ",s);
  for (i=0; i<numInputs+numOutputs; i++) {
    printf("%g ",p->data[i]);
  }
  printf("\n");
  return 0;
}


/**********************************************************************
 * Read all patterns from filename, place at end of Data linked list,
 * and reset data front and rear pointers.
 **********************************************************************/
void readPatterns (DataFile *dataFile) {
  FILE *file;
  int i,j, ret;
  Pattern *p;
  float *xp, junkFloat;
  //int greater=0, less=0;

  if((file = fopen(dataFile->name,"r")) == NULL) {
    printf("Couldn't open data file %s.\n",dataFile->name);
    exit(1);
  }

  for (i=0; i<dataFile->num; i++) {
    /* For each new pattern, add it to linked list of all patterns,
     * and save if first one in dataFile pointer */
    p = (Pattern *) malloc(sizeof(Pattern));
    p->data = (float *)malloc((numInputs+numOutputs)*sizeof(float)); /////
    if (i == 0) dataFile->first = p;
    if (firstPattern == NULL) 
      firstPattern = p;
    else 
      lastPattern->next = p;
    lastPattern = p;
    xp = p->data;
    for (j=0; j<numInputs+numOutputs; j++) {
      ret = fscanf(file, "%f", xp);
      if (ret == EOF) {
	printf("Premature end of file when reading %s.\n",dataFile->name);
	exit(1);
      }
      xp++;
    }
  }

  //fclose(file);

/**
  printf("File %s\n",dataFile->name);
  p = dataFile->first;
  for (i=0; i<dataFile->num; i++) {
    for (j = 0; j<numInputs+numOutputs; j++)
      printf("%g ",p->data[j]);
    printf("\n");
    p = p->next;
  }
**/

  if (fscanf(file,"%f",&junkFloat) != EOF) {
    printf("Warning: Unread data left in %s\n",dataFile->name);
  }

  fclose(file);
/*  fprintf(stderr,"Read %d lines from %s\n",dataFile->num,dataFile->name);*/
}

/**********************************************************************/

float *sumOutVars;
float *sumTargets;
float *sumSqTargets;

//void calcError(Set *set, float *wh1, float *wh2, float *wo,
//       float *rmsError, float *fracCorr, float *Rsquared, int print) {
void dumpOutput(Set *set, char *outputfile) {
  FILE *fout;
  int pattern, i, numCorr, numPatCorr;
  float *xp, *tp;

  if ((fout=fopen(outputfile, "w"))==NULL) {
    printf("Couldn't write to %s\n", outputfile);
    exit(1);
  }

  numCorr = 0;
  numPatCorr = 0;
  numwindows=0;
  for (pattern=0; pattern < set->num; pattern++) {

    xp = set->patterns[pattern]->data;
    tp = set->patterns[pattern]->data+numInputs;
    calcOutput(xp,wh1,wh2,wo);

    fprintf(fout, "pat %5d",pattern+1);
    fprintf(fout, " outputs");
    for (i=0; i<numOutputs; i++)  {
      fprintf(fout, " %10f",yo[i]);
      prob[numwindows][i]=yo[i];
    }
    fprintf(fout, "\n");
    numwindows++;
  } /* end of pattern loop */

  fclose(fout);
}
/**********************************************************************/

void calcOutput(float *x, float *wh1, float *wh2, float *wo) {
  int u, i; //, nio;
  float sum, *wp;

  for (u=0, wp=wh1; u<numHidden1; u++) {
    sum = 0.;
    for (i=0; i<numInputs; i++)
      sum += *wp++ * x[i];
    sum += *wp++;		/* bias weight */
    yh1[u] = logistic(sum);
  }

  if (numHidden2 > 0) {
    for (u=0, wp=wh2; u<numHidden2; u++) {
      sum = 0.;
      for (i=0; i<numHidden1; i++)
	sum += *wp++ * yh1[i];
      sum += *wp++;		/* bias weight */
      yh2[u] = logistic(sum);
    }
  }

  for (u=0, wp=wo; u<numOutputs; u++) {
    sum = 0.;
    if (numHidden2 > 0) {
      for (i=0; i<numHidden2; i++)
	sum += *wp++ * yh2[i];
    } else {
      if (numHidden1 > 0) {
	for (i=0; i<numHidden1; i++)
	  sum += *wp++ * yh1[i];
      }
    }
    if (numHidden1 == 0 || InputsThrough) {
      for (i=0; i<numInputs; i++)
	sum += *wp++ * x[i];
    }
    sum += *wp++;		/* bias weight */
    if (linearOut)
      yo[u] = sum;
    else
      yo[u] = logistic(sum);  /*sum;*/
  }
}

/**********************************************************************/

void
printNet(int ni, int nh1, int nh2, int no, float *x, float *wh1, float *wh2,
	 float *wo, float *yh1, float *yh2, float *yo)
{
  int u, i, nio, nh;

  if (nh2 > 0) 
    printf("Warning: Only printing first hidden layer and output layer.\n");

  if (nh1 > 0) {
 
    /* Printing first hidden layer */

    printf("\n    ");
    for (i=0; i<ni; i++) {
      printf("\n%4.2f",x[i]);
      for (u=0; u<nh1; u++)
	printf(" %6.2f",wh1[u*(ni+1)+i]);
    }
    printf("\n    ");
    for (u=0; u<nh1; u++)
      printf(" %6.2f",wh1[u*(ni+1)+ni]);

    printf("\n    ");
    for (u=0; u<nh1; u++)
      printf(" %6.2f",yh1[u]);

    /* Print second hidden layer */

    if (nh2 > 0) {

      printf("\n\n    ");
      nio = nh1+1;
      for (u=0; u<nh2; u++) {
	for (i=0; i<nh1; i++)
	  printf(" %6.2f",wh2[u*nh1+i]);
	printf(" %6.2f  %6.2f\n",wh2[u*nh1+nh1-1],yh2[u]);
      }
    }

    /* Print output layer */

    printf("\n\n    ");
    if (nh2 > 0) {
      nio = nh2 + 1;
      nh = nh2;
    } else {
      nio = nh1+1;
      nh = nh1;
    }
    if (InputsThrough)
      nio += ni;
    for (u=0; u<no; u++) {
      for (i=0; i< nh; i++)
	printf(" %6.2f",wo[u*nio+i]);
      if (InputsThrough) 
	for (i=0; i<ni; i++)
	  printf(" %6.2f",wo[u*nio+i+nh]);
      printf(" %6.2f  %6.2f\n",wo[u*nio+nio-1],yo[u]);
    }

  } else {

    /* just single output layer */

    printf("\n    ");
    for (i=0; i<ni; i++) {
      printf("\n%4.2f",x[i]);
      for (u=0; u<no; u++)
	printf(" %6.2f",wo[u*(ni+1)+i]);
    }
    printf("\n    ");
    for (u=0; u<no; u++)
      printf(" %6.2f",wo[u*(ni+1)+ni]);

    printf("\n    ");
    for (u=0; u<no; u++)
      printf(" %6.2f",yo[u]);
    printf("\n");
  }
}

/**********************************************************************/

float logistic(float x) {
  return (1. / (1. + exp(-x)));
}

/**********************************************************************/

float Dlogistic(float y) {
    return (y * (1 - y));
}

/**********************************************************************/
/* Return the number of rows and words, as counted by wc. */

void nLinesWords (char *filename, int *rows, int *words) {
  char string[100], temp_fname[30];
  FILE *f;
  int i;

  //*rows=4; *words=12; return;

  sprintf(temp_fname, "wcout_%s", filename);
  sprintf(string,"wc %s > %s",filename, temp_fname);
  i = mysystem(string);
  sleep(2);
  f = fopen(temp_fname,"r");
  fscanf(f,"%d %d",rows,words);
  fclose(f);
  sprintf(string, "rm -f %s", temp_fname);
  i=mysystem(string);
  sleep(2);
}

#include <sys/types.h>
#include <sys/time.h>

/*
		Addititive congruential  random number generator
		adapted from  Hill's Computer Graphics. Tests for 
		independence of sequence (Big-R) and uniformity of 
		distribution (Chi-square) are output. The histogram 
		dist. of random numbers is also output (3/22/91). 
		A floating point number between 0 - 1 is returned.
*/
 
/* global quantities for random number generator */
#define M 714025
#define IA 1366
#define IC 150889

long   the_seed;

/*--------------------------------------------------------------*/

void
init_acgrnd_by_time()
{
  time_t time();

  the_seed = time(NULL);
}

/*-------------------------------------------------------------*/

float acgrnd() {
      static long iy,ir[98];
      static int iff=0;
      int j;
      void nrerror();

      if (the_seed < 0 || iff == 0) {
              iff=1;
              if ((the_seed=(IC-(the_seed)) % M) < 0) the_seed = -(the_seed);
              for (j=1;j<=97;j++) {
                      the_seed=(IA*(the_seed)+IC) % M;
                      ir[j]=(the_seed);
              }
              the_seed=(IA*(the_seed)+IC) % M;
              iy=(the_seed);
    }
      j=1 + 97.0*iy/M;
      if (j > 97 || j < 1) exit(0);
      iy=ir[j];
      the_seed=(IA*(the_seed)+IC) % M;
      ir[j]=(the_seed);
      return (float) iy/M;
}

/**********************************************************************
  date() returns string like    Wed Jun 26 10:55:15 1991
**********************************************************************/

char *date()
{
  time_t time(),a;
  char *ctime(), *sfront; //, *send;

  a = time(NULL);
  sfront = ctime(&a);
  *(sfront+24) = '\0';
/*  send = sfront;
  for (send = sfront; *send != ' ' || *(send+1) != ' '; send++) ;
  *send = '\0';
*/
  return(sfront);
}
/**********************************************************************
  mysystem("command") does unix command and then sleeps for n seconds.
  returns exit status
**********************************************************************/
int mysystem(char *s) {
  int i;
  if (debug) {printf("Doing command %s\n",s); fflush(stdout);}
  i = system(s);
/*  sleep(5);*/
  return i;
}
/**********************************************************************
  fexists("file.txt")   returns 1 if "file.txt" can be opened 'r', 0 otherwise
**********************************************************************/

int fexists(char *filename) {
  FILE *fid;

  if (debug) {printf("Does %s exist?",filename); fflush(stdout);}

  if ((fid=fopen(filename,"r")) != NULL) {
    fclose(fid);
    if (debug) { printf("  yes.\n"); fflush(stdout);}
    return 1;
  } else {
    if (debug) { printf("  no.\n"); fflush(stdout); }
    return 0;
  }
}


