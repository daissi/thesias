/*
 * This file is part of THESIAS
 * Copyright (C) 2004-2020 David-Alexandre Trégouët, Valérie Garelle
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// LA FONCTION INVERSION A ETE REMPLACEE PAR SYSL ET LES DOUBLE POINTEURS PAR DES TABLEAUX DE DIM 2


/* Tout se programme sera réalisé pour des individus indépendants */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>

//#include "mconf.h"
//#include "definition.h"


//DEB debut de definition.h
/*********************************************************************************************/
/*						      Déclarations des Constantes propres à SEM  					 */
/***************************************************************************************** ****/


#define  mymaxit  2501/*2501*/    /*nb d'itération*/
#define  nburn  1500/*1500*/     /*nb de burns*/
#define  seuilfisher 0.0000001 /*seuil de convergence de l'algorithme de fisher*/
#define  maxrun 20      /*nombre d'iteration de l'aglorithme de ' */

#define maxn		100  /*maximum nb of estimated parameters*/
#define pi 3.141592654
#define  maxloc 25  	/*nb maximum de loci */
#define  locmq  2       /*nb de genotype manquant autorisé < nbloci-locmq*/
#define  maxcov 12       /*nb maximum de covariables and interaction*/
#define  maxhypor 4     /*maximal number of different hypotheses on OR/differences */
#define  maxor   8      /*maximal number of OR/differences to be tested within the same hypothesis*/
#define  maxconst 15    /*maximal number of constraints to be set */
#define  maxadd 7       /*maximal number of additivity constraints to be set */
#define  maxhypint 6
#define maxclass 4

#define TINY 1.0e-20;
// Fin definition.h




extern double igam ( double, double );
extern double igami ( double, double );
extern double polevl ( double, double *, int );
extern double ndtri (double);
extern double gamma (double);
extern double chdtrc(double, double);
extern int mtherr(char *, int );
// VG
//extern void affichage(FILE *,double , double ),affichage2(FILE *,double , double ),affich3(FILE *,double , double );
//extern void fisherscoring(), Xfisherscoring();
//extern void fishem(double *, double *, matrixp ),Xfishem(double *, double *, matrixp );
//extern double likelihood(double *, double *) , Xlikelihood(double *, double *);
//extern double condlike(double *), Xcondlike(double *);
//extern void phenomean(FILE *,FILE *,matrixp);
//extern void tritime(),tripair(),coxtempo();
//extern void initlist();
//extern void coxrubin(double **);
//extern void likematchpair(double *,double *,double *);
//extern void breslow1(double *,double *,double *);
//extern void matchpair();
//extern void fishpair(double *, double *,matrixp);
//extern int ipow(int ,int);
//extern double residuel(double *, double *);
//extern void presence(),hapopres(),recodage(),distrmq();
//extern double sysl(double **, long );
//extern int coding(double);
//extern void rsquare(double *,FILE *,FILE *);
//extern void Xfishnull(double *,matrixp),fishnull(double *,matrixp);
//extern void ludcmp(double **,int, int *, double *),lubksb(double **,int,int *, double *), inversion(double **,int );
//FIN VG

//VG
//extern void pairedfish(double **);
// FIN VG

// DEB definition.h


/*********************************************************************************************/
/*						      Déclarations des Types    									 */
/*********************************************************************************************/

typedef int  vectgen[maxloc][2];
typedef double vectdon[maxcov];

typedef struct dhaplotype { /* Definition d'un haplotype */
	int numnew,identif;             /* Verifier comment est utilise identif */
	short present, listall[maxloc]; /* Verifier si "present" est utilisé */
	char code[maxloc][2];           /* Verifier si code est utilisé qu'à la fin.
	Si oui,on pourrait ne pas le définir en global.*/
	double frqle;
	struct dhaplotype *down;
} dhaplotype;


typedef struct dindividu {  /* Definition d'un individu */
	int tnbhapo,nblm; /* nb de paires haplotypiques possibles + nb de loci manquants */
    int **idnb;/* ensemble des paires haplotypiques possibles*/
    int hapest[2]; /* paire haplotypique attribuee */
	vectgen marq; /* vecteur des genotypes*/
	double phen[2],wgt; /*vecteur phenotypique */
	vectdon z; /* vecteur des covariables */
	//int gender;
	struct dindividu *next;
} dindividu;


typedef struct combgeno
{double fquence;
 int settab[maxloc],place;
 int **sethap,nbcomb;
 struct combgeno *apres;
} combgeno;



typedef struct dliste {
	double *compo;
	double pival;
	struct dliste *proch;
 } dliste;



/*FISHER SCORING*/
typedef double vectp[maxn];
typedef vectp matrixp[maxn];


/*********************************************************************************************/
/*						      Déclarations des variables 									 */
/*********************************************************************************************/

double **tablo,*tabpi,*dmat2,*mdvs2,*modif2,mdvd2[maxn][maxn];
double *tabres;

//dliste *listbase,*listp,*listold;

dhaplotype *tnbhnew,*tnbhbase,*vect1,*vect2;
dindividu *base,*suiv;

int nbtotused;
int nbhhypo; /* nb potentiel d'haplotypes*/
int nbloci;  /* nb de loci étudies */
int maxhapair=0;/* maximum nb de paires haplotypiques possibles observees */
int nnt;     /* nb of regression parameters  to be estimated */
int *nbor;
int ***tabinter,**tabhypo,**tadd,**tabint,**tabhypint;

short *fcoda1; /*Indicateur d' eventuelle presence (1/0) des haplotypes */
int  *fcoda2; /*Tableau de renumerotation des haplotypes eventuellement presents */

short msdata;  /* indicateur de prise en compte des donnees manquantes*/
short ldmatrix;
short letcod;  /* indicateur de prise en compte de code pour les alleles*/
short chxt,offset;    /* type of phenotype: binary or quantitative or survival */
int ajust;   /* nb de variables d'ajustement */
short rsq,wald;
int numajust[maxcov],idoffset;
short nbcas,nbcasm,nbtem,nbtot,nbused; /* nb de cas avec et sans missing data, total de sujets */
short nbhf[2][3];
short haplozero;  /* indicateur d'estimation des effets haplotypiques */
//int interor,hypoth,nbhypor,nbadd=0,hypint=0,intercov=0,xlnk=0;
//VG 09112006
int interor=0,hypoth=0,nbhypor=0,nbadd=0,hypint=0,intercov=0,xlnk=0;

double *alfreq,*hafreq; /*tableau des fqces alléliques et haplotypiques sous independances*/
double *freqest, *tempfreq; /*tableau des fqces haplotypiques estimées et courantes*/
double **freqdist, **tempdist; /*tableau des fqces haplotypiques estimées et courantes chez les temoins/cas*/
int *tabmq;     /* tableau des donnees manquantes */
char letter[maxloc][2]; /* tableau du codage des alleles */
double *moyeff,*effest,**vecbeta,**vecbetat,*moyfreq;

double tempx;
short *inclus;
int *numhap,nbhest=0;
double mean,ste,ste0,steres,meanste=0;
double llambdaval,wbgama,wblamda;


int nbcatego,nkat;
int *nbsujktgo;

FILE *parafile;

enum {go, quit, restart} conti;

int *itp,*nitp,*itptp,*nitptp;
int nall, n;

//FIN definition.h


/*********************************************************************************************/
/*						      Déclarations des fonctions 									 */
/*********************************************************************************************/
//void lecture();
void determhapo();
void nbhapo1(vectgen);
void nbhapo0(vectgen);
void initfreq();
void lecteffe();
void generhap();
double probatot();
double probacond(int);
double Xprobatot();
double Xprobacond(int);

double llambda(double );

double somdelai();


void polytomous();
void categorie();
void vpolyto(double **);
double likepoly(double *,double *);


//deb foncgener

void ludcmp(double **a,int n,int *indx,double *d)
{
 int  i, imax,j,k;
 double big,dum,sum,temp;
 double *vv;
 vv=(double *)malloc((size_t) (n*sizeof(double)));
 *d=1.0;
 for (i=0;i<n;i++)
 {big=0.0;
  for (j=0;j<n;j++) if ((temp=fabs(a[i][j]))>big) big=temp;
  if (big==0.0) printf("Singular matrix in routine LUDCMP\n");
  vv[i]=1.0/big;
 }
 for (j=0;j<n;j++)
 {for (i=0;i<j;i++)
  {sum=a[i][j];
   for (k=0;k<i;k++) sum-=a[i][k]*a[k][j];
   a[i][j]=sum;
  }
  big=0.0;
  for (i=j;i<n;i++)
  {sum=a[i][j];
   for (k=0;k<j;k++) sum-=a[i][k]*a[k][j];
   a[i][j]=sum;
   if ((dum=vv[i]*fabs(sum))>=big) {big=dum;imax=i;}
  }
  if (j!=imax)
  {for (k=0;k<n;k++)
   {dum=a[imax][k];
    a[imax][k]=a[j][k];
	a[j][k]=dum;
   }
   *d=-(*d);
   vv[imax]=vv[j];
  }
  indx[j]=imax;
  if (a[j][j]==0.0) a[j][j]=TINY;
  if (j!=n-1)
  {dum=1.0/(a[j][j]);
   for (i=j+1;i<n;i++) a[i][j] *=dum;
  }
 }
 vv=NULL; free (vv);
}



void lubksb(double **a,int n,int *indx,double *b)
{int i,ii=0,ip,j;
 double sum;
 for (i=1;i<=n;i++)
 {ip=indx[i-1];
  sum=b[ip];
  b[ip]=b[i-1];
  if (ii) for (j=ii;j<=i-1;j++) sum-=a[i-1][j-1]*b[j-1];
  else if (sum) ii=i;
  b[i-1]=sum;
 }
 for (i=n;i>=1;i--)
 {sum=b[i-1];
  for (j=i+1;j<=n;j++) sum-=a[i-1][j-1]*b[j-1];
  b[i-1]=sum/a[i-1][i-1];
 }
}

void inversion(double **a,int n)
{int *indx,i,j;
 double d,*col,**y;
 indx=(int *) malloc((size_t) (n*sizeof(int)));
 col=(double *) malloc((size_t) (n*sizeof(double )));
 y=(double **) malloc((size_t) (n*sizeof(double *)));
 for (i=0;i<n;i++) y[i]=(double *) malloc((size_t) (n*sizeof(double)));


  ludcmp(a,n,indx,&d);

 for (j=0;j<n;j++)
 {for (i=0;i<n;i++) col[i]=0.0;
  col[j]=1.0;
  lubksb(a,n,indx,col);
  for (i=0;i<n;i++) y[i][j]=col[i];
 }
 for (i=0;i<n;i++) for (j=0;j<n;j++) a[i][j]=y[i][j];

 indx=NULL;free (indx);
 col=NULL;free (col);
 y=NULL;free (y);

}

/*
LE COMMENTER UNIQUEMENT POUR WIN
*/

int ipow(int val,int expv)
 {int i,res;
  res=1;
  for (i=0;i<expv;i++) res*=val;
  return(res);
 }
/*FIN WIN!*/


/********************Identifie les haplotypes potentiels du fichiers*************************/
void hapopres()
{vect1=tnbhbase;
 nbhhypo=0;
 while (vect1!=NULL)
 {vect1->present=0;vect1->identif=-1;
  if (fcoda1[vect1->numnew]==1)
		  {vect1->present=1;vect1->identif=nbhhypo;
	       fcoda2[vect1->numnew]=nbhhypo;
		   nbhhypo+=1;
		  }
		  vect1=vect1->down;
	 }
 vect1=NULL;

}

/*******************Renumerote les haplotypes des individus**********************************/
void recodage()
{int i;
 suiv=base;
 do
 {for (i=0;i<suiv->tnbhapo;i++)
  {suiv->idnb[i][0]=fcoda2[suiv->idnb[i][0]];
   suiv->idnb[i][1]=fcoda2[suiv->idnb[i][1]];
  }
  suiv=suiv->next;
 }while ((suiv!=NULL) && (suiv->next!=NULL));
 suiv= NULL;
}

/**********************************Calculs du nombre de donnees manquantes********************/
void distrmq()
{int i,j;
 double stdd=0;
 for (i=0;i<2;i++) for (j=0;j<3;j++) nbhf[i][j]=0;
 mean=0;
 nbused=0;
 nbcas=0;nbtot=0;nbcasm=0;
 tabmq=(int *) malloc((size_t) ((nbloci+1)*sizeof(int)));
 for (i=0;i<nbloci+1;i++) tabmq[i]=0;
 suiv=base;
 if (xlnk==0)
 {do
  {tabmq[suiv->nblm]+=1;
   if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
   {nbcas+=(suiv->phen[0]==1)*((msdata==0)*(suiv->nblm==0)+(msdata==1)*(suiv->nblm<nbloci-locmq+1));
    nbcasm+=(suiv->phen[0]==1)*(suiv->nblm==0);
   }
   if (chxt==2)
   {mean+=(suiv->phen[0])*((msdata==1)*(suiv->nblm<nbloci-locmq+1)+(msdata==0)*(suiv->nblm==0));
    stdd+=(suiv->phen[0])*(suiv->phen[0])*((msdata==1)*(suiv->nblm<nbloci-locmq+1)+(msdata==0)*(suiv->nblm==0));
   }
   nbused+=1*((msdata==1)*(suiv->nblm<nbloci-locmq+1)+(msdata==0)*(suiv->nblm==0));
   nbtot+=1;
   suiv=suiv->next;
  } while ((suiv!=NULL) && (suiv->next!=NULL));
 }
 else if (xlnk==1)
 {do
  {if (suiv->nblm==0)   nbhf[(int) suiv->z[0]][0]++;
   tabmq[suiv->nblm]+=1;
   if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
   {nbcas+=(suiv->phen[0]==1)*((msdata==0)*(suiv->nblm==0)+(msdata==1)*(suiv->nblm<nbloci-locmq+1));
    nbcasm+=(suiv->phen[0]==1)*(suiv->nblm==0);
    /* ATTENTION QUAND IL Y AURA DES DONNES MANQUANTES AVEC CETTE OPTION*/
    if (suiv->nblm==0) nbhf[(int) suiv->z[0]][1+(int) suiv->phen[0]]++;
   }
   if (chxt==2)
   {mean+=(suiv->phen[0])*((msdata==1)*(suiv->nblm<nbloci-locmq+1)+(msdata==0)*(suiv->nblm==0));
    stdd+=(suiv->phen[0])*(suiv->phen[0])*((msdata==1)*(suiv->nblm<nbloci-locmq+1)+(msdata==0)*(suiv->nblm==0));
   }
   nbused+=1*((msdata==1)*(suiv->nblm<nbloci-locmq+1)+(msdata==0)*(suiv->nblm==0));
   nbtot+=1;
   suiv=suiv->next;
  }while ((suiv!=NULL) && (suiv->next!=NULL));
 }

 suiv=NULL;
 if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6)) {nbtem=nbused-nbcas;}
 ste=sqrt((stdd-mean*mean/nbused)/(nbused-1));
 ste0=ste;
 mean/=nbused;
 }


 int coding(double val)
{int i,j;
 short trouve=0;
 j=0;i=-1;
 while ((j<nbhhypo) && (trouve==0))
 {if (val==numhap[j])
  {trouve=1;i=j;}
  j++;
 }
 if (trouve==0) {i=-1;}
 return(i);
}


void affichage(FILE *outf,double valef, double valse)
{double valt, valp;
 valt=(valef/valse);
 valp=valt*valt;
 fprintf(outf,"<td align=left>%f</td><td align=left>%f</td></tr>\n",valse,valt);
 if ( (chxt==1) || (chxt==4))
 {fprintf(outf,"<tr><td align=right colspan=""5"">OR = %.5f [%.5f - %.5f] ",exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));
 }
 else if (chxt==2)
 {fprintf(outf,"<tr><td align=right colspan=""5"">Diff = %.5f [%.5f - %.5f] ",valef,valef-1.96*valse,valef+1.96*valse);
 }
 else if ((chxt==3) || (chxt==5))
 {fprintf(outf,"<tr><td align=right colspan=""5"">HRR = %.5f [%.5f - %.5f] ",exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));
 }
 if (valp>0)
 {fprintf(outf,"p=%f </td></tr>\n\n",chdtrc(1.0,valp));
 }
 else fprintf(outf,"p is undefined </td></tr>\n\n");
}

void affichage2(FILE *outf,double valef,double valse)
{double valt, valp;
 valt=(valef/valse);
 valp=valt*valt;
 fprintf(outf,"%f\t%f\n",valse,valt);
 if ( (chxt==1) || (chxt==4))
 {
	 //fprintf(outf,"\t\t\tOR = %.5f [%.5f - %.5f] ",exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));
	 fprintf(outf,"\t\tOR = %.5f [%.5f - %.5f] ",exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));//<me>
 }
 else if (chxt==2)
 {
	 //fprintf(outf,"\t\t\tDiff = %.5f [%.5f - %.5f] ",valef,valef-1.96*valse,valef+1.96*valse);
	 fprintf(outf,"\t\tDiff = %.5f [%.5f - %.5f] ",valef,valef-1.96*valse,valef+1.96*valse);//<me>
 }
 else if ((chxt==3) || (chxt==5))
 {
	 //fprintf(outf,"\t\t\tHRR = %.5f [%.5f - %.5f] ",exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));
	 fprintf(outf,"\t\tHRR = %.5f [%.5f - %.5f] ",exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));//<me>
 }
 if (valp>0)
 {fprintf(outf," p=%f\n\n",chdtrc(1.0,valp));
 }
 else fprintf(outf,"  p is undefined \n\n");
}

void affich3(FILE *outf,double valef, double valse)
{double valt, valp;
 valt=(valef/valse);
 valp=valt*valt;



 if ((valef==0) && ( (valse==1) || (valse==0)))
 {if ( (chxt==1) || (chxt==4)) {fprintf(outf,"OR = 1 ");}
  else if (chxt==2) {fprintf(outf,"Diff = 0 ");}
 }
 else if (valse>0)
 {if ( (chxt==1) || (chxt==4)) {fprintf(outf,"OR = %.5f [%.5f - %.5f] ", exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));}
  else if (chxt==2)
  {fprintf(outf,"Diff = %.5f [%.5f - %.5f] ",valef,valef-1.96*valse,valef+1.96*valse);}
  else if ((chxt==3) || (chxt==5)) {fprintf(outf,"HRR = %.5f [%.5f - %.5f] ", exp(valef),exp(valef-1.96*valse),exp(valef+1.96*valse));}
  if (valp>0) {fprintf(outf," p=%f\n",chdtrc(1.0,valp));}
  else fprintf(outf,"  p is undefined \n");
 }
 else if (valef!=0)
 {if ( (chxt==1) || (chxt==4)) {fprintf(outf,"OR = %.5f ", exp(valef));}
  else if (chxt==2)
  {fprintf(outf,"Diff = %.5f ",valef);
  }
 }
}

void presence()
{int i;
 suiv = base;
 for (i=0;i<nbhhypo;i++) inclus[i]=0;
  do
 {if (suiv->tnbhapo>0)
  {inclus[suiv->hapest[0]]=1;
   inclus[suiv->hapest[1]]=1;
  }
  suiv = suiv->next;
 } while ((suiv!=NULL) && (suiv->next!=NULL));

}


double sysl(double a[maxn][maxn], long nn)
{/* Inversion et calcul du determinant d'une matrice */
 long ia[maxn];
 long i, j,kk, ii;
 double cf, d,val;
 int su1, su2;

  /*sysl*/
  val = 1;
  kk = 1;
  while (kk <= nn)
  {ia[kk - 1] = 0;
   i = 1;
   su1 = 1;
   while ((i <= nn) && (su1==1))
   {j = 1;
    su2 = 1;
    while ((j <= kk) && (su2==1))
	{if (i == ia[j - 1]) {su2 = 0;}
	 j++;
    }
    ii = i;
    i++;
    if (su2==1)
	{if (a[ii - 1][ii - 1] != 0.0)
     {ia[kk - 1] = ii;
	  d = a[ii - 1][ii - 1];
	  su1 = 0;
	 }
    }
   }
   if (su1==1)
   {printf("Matrice non inversible\n");
    val=0;
	exit(0);
   }
   for (i = 0; i < nn; i++)
   {if (i + 1 != ii)
    {cf = a[i][ii - 1] / d;
	 for (j = 0; j < nn; j++)
	 {if (j + 1 != ii) a[i][j] -= a[ii - 1][j] * cf;
	 }
	}
   }
   for (j = 0; j < nn; j++)
   {if (j + 1 != ii)
    {a[ii - 1][j] /= d;
     a[j][ii - 1] = -(a[j][ii - 1] / d);
	}
   }
   a[ii - 1][ii - 1] = 1.0 / d;
   val *= d;
   kk++;
  }
  return(val);
}



void frqcomb(combgeno *genotb,double *moyfreq)
{int i,nbhet,hetero[maxloc],h0,h00,tnbh1,backha,numh,idp,idx;
 /*int **listhap;*/ short iref;
 double answ;
 tnbh1=0;
 nbhet=0;
 for (i=0;i<nbloci;i++)
 {if (genotb->settab[i]==1)
  {hetero[nbhet]=i+1;nbhet+=1;}
 }
 if (nbhet<2)
 {tnbh1=1;
  h0=0;h00=0;
  for (i=1;i<=nbloci;i++)
  {h0+=(genotb->settab[i-1]==2)*ipow(2,nbloci-i);
   h00+=(genotb->settab[i-1]>0)*ipow(2,nbloci-i);
  }
  h0=fcoda2[h0];
  h00=fcoda2[h00];
  answ=moyfreq[h0]*moyfreq[h00]*(2-(h0==h00));
  genotb->sethap= (int **) malloc((size_t) (1*sizeof(int *)));
  for (i=0;i<tnbh1;i++) genotb->sethap[i]=(int *)malloc((size_t) (2*sizeof(int)));
  genotb->sethap[0][0]=h0;
  genotb->sethap[0][1]=h00;
 }
 else
 {tnbh1=ipow(2,nbhet-1);
  genotb->sethap= (int **) malloc((size_t) (tnbh1*sizeof(int *)));
  for (i=0;i<tnbh1;i++) genotb->sethap[i]=(int *)malloc((size_t) (2*sizeof(int)));
  //listhap=(int **)malloc((size_t) (tnbh1*sizeof(int *)));
  //for (i=0;i<tnbh1;i++) listhap[i]=(int *)malloc((size_t) (2*sizeof(int)));
  backha=0;
  for (i=1;i<=nbloci;i++)
  {if (genotb->settab[i-1]!=1) {backha+=(genotb->settab[i-1]==2)*ipow(2,nbloci-i);}
  }
  for (i=0;i<tnbh1;i++)
  {genotb->sethap[i][0]=backha;
   genotb->sethap[i][1]=backha+ipow(2,nbloci-hetero[0]);
  }

  for (i=1;i<=nbhet-1;i++)
  {numh=0;idx=ipow(2,nbhet-1-i);iref=1;idp=1;
   do
   {genotb->sethap[numh][0]+=(iref==2)*ipow(2,nbloci-hetero[i]);
    genotb->sethap[numh][1]+=(iref==1)*ipow(2,nbloci-hetero[i]);
    idp++;
    if (idp>idx) {if (iref==1) {iref=2;}
                  else         {iref=1;}
                  idp=1;
			     }
    numh++;
    }while (numh<tnbh1);
  }
  answ=0;
  for (i=0;i<tnbh1;i++)
  {genotb->sethap[i][0]=fcoda2[genotb->sethap[i][0]];
   genotb->sethap[i][1]=fcoda2[genotb->sethap[i][1]];
   h0=genotb->sethap[i][0];h00=genotb->sethap[i][1];
   answ+=moyfreq[h0]*moyfreq[h00]*(2-(h0==h00));
  }
 // listhap=NULL;
 // free (listhap);
 }
 genotb->nbcomb=tnbh1;
 genotb->fquence=answ;

}


void rsquare(double *moyfreq,FILE *outfile,FILE *outres)
{int ixpar,i,trouve,nctab=0,ll,j,h0,h00;
 combgeno *genotabase,*genotab1,*genotab2;
 double *expect,numera,denomi;
 char l;

 genotabase=(combgeno *)malloc((size_t) (sizeof(combgeno)));
/* ATTENTION AU DONNES MAQNAUTENS Pour l'instant supposant qu'il n'y en a pas*/
 ixpar=0;
 genotab1=genotabase;
 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {if (suiv->nblm==0)
  {if (ixpar==0)
   {for (i=0;i<nbloci;i++) genotab1->settab[i]=(suiv->marq[i][0]+suiv->marq[i][1]-2);
    genotab1->place=nctab;
    frqcomb(genotab1,moyfreq);
    genotab1->apres=(combgeno *)malloc((size_t) (sizeof(combgeno))); genotab1=genotab1->apres;
    ixpar=1; nctab++;
   }
   else
   {genotab2=genotabase;
    trouve=0;
    while ((genotab2!=genotab1) && (trouve==0))
    {i=0;while ((genotab2->settab[i]==(suiv->marq[i][0]+suiv->marq[i][1]-2)) && (i<nbloci))	{i++;}
	 if (i==nbloci) {trouve=1;} else trouve=0;
	 genotab2=genotab2->apres;
    }
    if (trouve==0)
    {for (i=0;i<nbloci;i++) genotab1->settab[i]=(suiv->marq[i][0]+suiv->marq[i][1]-2);
     genotab1->place=nctab;
     frqcomb(genotab1,moyfreq);
	 genotab1->apres=(combgeno *)malloc((size_t) (sizeof(combgeno)));
	 genotab1=genotab1->apres;
	 nctab++;
    }
   }
  }
  suiv = suiv->next;
 }
 genotab1->apres=NULL;
 genotab1=NULL;


 expect=(double *)malloc((size_t) (nbhhypo*sizeof(double)));
 for (i=0;i<nbhhypo;i++) expect[i]=0.0;

 fprintf(outfile,"<tr><td align=left width=20%%></td></tr>\n");
 fprintf(outfile,"<tr><td align=left width=20%%></td></tr>\n");
 fprintf(outfile,"<tr>\n<td align=center colspan=""5"">Haplotypic R2 information</td></tr>\n\n");
 fprintf(outfile,"<tr><td align=left width=20%%></td></tr>\n");
 fprintf(outfile,"<tr><td align=left> </td><td align=left></td><td align=left>Frequency</td><td align=left>R Square</td>\n");
 fprintf(outfile,"<td> </td></tr>\n\n");


 fprintf(outres,"\t\tHaplotypic R2 information\n\n");

 vect1=tnbhbase;i=0;
 while (vect1!=NULL)
 {if (vect1->present==1)
  {ll=fcoda2[vect1->numnew];
   if (moyfreq[ll]>0)
   {//printf("I=%d Fq=%lf\t",ll,moyfreq[ll]);
	expect[ll]=0.0;
    genotab1=genotabase;
    while (genotab1->apres!=NULL)
 	{numera=0.0;denomi=0.0;
     for (j=0;j<genotab1->nbcomb;j++)
	 {h0=genotab1->sethap[j][0];h00=genotab1->sethap[j][1];
	  denomi+=moyfreq[h0]*moyfreq[h00];
      numera+=moyfreq[h0]*moyfreq[h00]*((ll==h0)+(ll==h00));
	 }
     expect[ll]+=genotab1->fquence*(numera/denomi)*(numera/denomi);
	 genotab1=genotab1->apres;
	}
	genotab1=NULL;
    expect[ll]-=4*moyfreq[ll]*moyfreq[ll];
    expect[ll]/=2*moyfreq[ll]*(1-moyfreq[ll]);
	//printf("R2 =%lf\n",expect[ll]);

	fprintf(outfile,"<tr><td align=left width=20%%>Haplotype [%d] </td>\n",i);
    fprintf(outfile,"<td align=center width =25%%>");
    fprintf(outres,"Haplotype [%d] \t",i);
	for (j=0;j<nbloci;j++)
    {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
     fprintf(outfile,"%c",l);
     fprintf(outres,"%c",l);
    }
	fprintf(outfile,"</td><td align=left>");
    fprintf(outfile,"%f</td>",moyfreq[ll]);
	fprintf(outres," %f\t",moyfreq[ll]);
	fprintf(outfile,"<td align=left>R2 = %f</td><td align=left>&nbsp;</td></tr>\n",expect[ll]);
    fprintf(outres,"R2 = %f\t\n",expect[ll]);

   }
  }
  i++;
  vect1=vect1->down;
 }
  /* NEW DAVID */
  free (genotabase);genotabase=NULL;
  free (genotab1);genotab1=NULL;
  free (genotab2);genotab2=NULL;
  free (expect);expect=NULL;


}



void fishnull(double *frqsem,matrixp sfisher)
{int i,j,dimf,dimf2;
 double h1,h2,like,likeli,val3;
 double *derivfrq,*deriv;int *place;

 double fisher[maxn][maxn];

 printf("Running Variance Estimation\n");
 place=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
 for (i=0;i<nbhhypo;i++) place[i]=-1;
 j=0;for (i=0;i<nbhhypo;i++) if (frqsem[i]>0) {place[i]=j;j++;}


 dimf=j-1;dimf2=j;
 //fisher=(double **) malloc((size_t) (dimf*sizeof(double *)));
 //for (i=0;i<dimf;i++) fisher[i]=(double *) malloc((size_t) (dimf*sizeof(double)));
 for (i=0;i<dimf;i++) for (j=0;j<dimf;j++) fisher[i][j]=0;

 derivfrq=(double *) malloc((size_t) (dimf2*sizeof(double)));
 for (i=0;i<dimf2;i++) derivfrq[i]=0;

 deriv=(double *) malloc((size_t) (dimf*sizeof(double)));
 for (i=0;i<dimf;i++) deriv[i]=0;

 suiv=base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
 {for (i=0;i<dimf2;i++) derivfrq[i]=0;
  for (i=0;i<dimf;i++) deriv[i]=0;
  /*CALCUL DE P(Y,G)*/
  likeli=0;
  for (i=0;i<suiv->tnbhapo;i++)
  {like=0;
   h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
   if ( (h1>0) && (h2>0))
   {like=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
    /*Calcul des dérivees premieres par rapport aux frequences */
    if (suiv->idnb[i][0]==suiv->idnb[i][1])
    {derivfrq[place[suiv->idnb[i][0]]]+=2*frqsem[suiv->idnb[i][0]];}
    else
    {derivfrq[place[suiv->idnb[i][0]]]+=2*frqsem[suiv->idnb[i][1]];
     derivfrq[place[suiv->idnb[i][1]]]+=2*frqsem[suiv->idnb[i][0]];
    }
    likeli+=like;
   }/* if */
  }/* for i*/
  /* ATTENTION A LA PREMIERE FREQUENCE D'OU LE -1 */

  if (likeli>0)
  {for (i=0;i<dimf;i++) {deriv[i]=derivfrq[i+1]/likeli;}
   for (i=0;i<dimf;i++) for (j=0;j<dimf;j++) fisher[i][j]+=deriv[i]*deriv[j];
  }
  suiv=suiv->next;
 }/* while */


 printf("Inverting Variance Matrix....\n");
 sysl(fisher,dimf);
 for (i=0;i<dimf;i++) for (j=0;j<dimf;j++) sfisher[i][j]=fisher[i][j];

 deriv=NULL;derivfrq=NULL;place=NULL;
 //fisher=NULL;free (  fisher);
 free ( deriv);free (  derivfrq);free( place);
}



void Xfishnull(double *frqsem,matrixp sfisher)
{int i,j,dimf,dimf2;
 double h1,h2,like,likeli,val3;
 double *derivfrq,*deriv;int *place;

 double fisher[maxn][maxn];

 printf("Running Variance Estimation\n");
 place=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
 for (i=0;i<nbhhypo;i++) place[i]=-1;
 j=0;for (i=0;i<nbhhypo;i++) if (frqsem[i]>0) {place[i]=j;j++;}


 dimf=j-1;dimf2=j;
 //fisher=(double **) malloc((size_t) (dimf*sizeof(double *)));
 //for (i=0;i<dimf;i++) fisher[i]=(double *) malloc((size_t) (dimf*sizeof(double)));
 for (i=0;i<dimf;i++) for (j=0;j<dimf;j++) fisher[i][j]=0;

 derivfrq=(double *) malloc((size_t) (dimf2*sizeof(double)));
 for (i=0;i<dimf2;i++) derivfrq[i]=0;

 deriv=(double *) malloc((size_t) (dimf*sizeof(double)));
 for (i=0;i<dimf;i++) deriv[i]=0;

 suiv=base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {for (i=0;i<dimf2;i++) derivfrq[i]=0;
  for (i=0;i<dimf;i++) deriv[i]=0;
  /*CALCUL DE P(Y,G)*/
  likeli=0;
  if ((int) suiv->z[0]==1)
  {for (i=0;i<suiv->tnbhapo;i++)
   {like=0;
    h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
    if ( (h1>0) && (h2>0))
    {like=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
     /*Calcul des dérivees premieres par rapport aux frequences */
     if (suiv->idnb[i][0]==suiv->idnb[i][1])
     {derivfrq[place[suiv->idnb[i][0]]]+=2*frqsem[suiv->idnb[i][0]];}
     else
     {derivfrq[place[suiv->idnb[i][0]]]+=2*frqsem[suiv->idnb[i][1]];
      derivfrq[place[suiv->idnb[i][1]]]+=2*frqsem[suiv->idnb[i][0]];
     }
     likeli+=like;
    } /* if */
   }/* for i*/
   /* ATTENTION A LA PREMIERE FREQUENCE D'OU LE -1 */
  }
  else if ((int) suiv->z[0]==0)
  {for (i=0;i<suiv->tnbhapo;i++)
   {like=0;
    h1=frqsem[suiv->idnb[i][0]];
    if (h1>0)
    {like=h1;
     /*Calcul des dérivees premieres par rapport aux frequences */
     derivfrq[place[suiv->idnb[i][0]]]+=1;
     likeli+=like;
    } /* if */
   }
  }
  if (likeli>0)
  {for (i=0;i<dimf;i++) {deriv[i]=derivfrq[i+1]/likeli;}
   for (i=0;i<dimf;i++) for (j=0;j<dimf;j++) fisher[i][j]+=deriv[i]*deriv[j];
  }
  suiv=suiv->next;
 }/* while */


 printf("Inverting Variance Matrix....\n");
 sysl(fisher,dimf);
 for (i=0;i<dimf;i++) for (j=0;j<dimf;j++) sfisher[i][j]=fisher[i][j];

 deriv=NULL;derivfrq=NULL;place=NULL;
 //fisher=NULL;free (  fisher);
 free ( deriv);free (  derivfrq);free( place);
}

//fin foncgener


//deb fonc12.c


void tritime()
{dindividu *parcour,*newsui,*proch;
 int trouve;

 suiv=base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {
  trouve=0;
  parcour=base;
  while ( (parcour!=NULL) && (parcour!=suiv) && (trouve==0))
  {if (suiv->phen[1]>=parcour->phen[1])
   {trouve=1;break;}
   parcour=parcour->next;
  }
  if (trouve==1)
  {newsui=base;
   proch=suiv->next;
   if (parcour==base)
   {newsui=base;
    while (newsui->next!=suiv) {newsui=newsui->next;}
    newsui->next=proch;
	base=suiv;
    suiv->next=parcour;
   }
   else
   {newsui=base;
    while (newsui->next!=parcour)
	{newsui=newsui->next;}
     newsui->next=suiv;

	newsui=parcour;
	while (newsui->next!=suiv)
	{newsui=newsui->next;}
     newsui->next=proch;

	suiv->next=parcour;
   }
   suiv=proch;
  }
  else suiv=suiv->next;
 }
 suiv= NULL;
proch= NULL; parcour=NULL; newsui=NULL;
free (newsui);free (proch);free(parcour);
}


void tripair()
{dindividu *parcour,*newsui,*proch;
 int trouve;

 suiv=base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {
  trouve=0;
  parcour=base;
  while ( (parcour!=NULL) && (parcour!=suiv) && (trouve==0))
  {if (suiv->phen[1]==parcour->phen[1])
   {trouve=1;break;}
   parcour=parcour->next;
  }
  if (trouve==1)
  {newsui=base;
   proch=suiv->next;
   if (parcour==base)
   {newsui=base;
    while (newsui->next!=suiv) {newsui=newsui->next;}
    newsui->next=proch;
	base=suiv;
    suiv->next=parcour;
   }
   else
   {newsui=base;
    while (newsui->next!=parcour)
	{newsui=newsui->next;}
     newsui->next=suiv;

	newsui=parcour;
	while (newsui->next!=suiv)
	{newsui=newsui->next;}
     newsui->next=proch;

	suiv->next=parcour;
   }
   suiv=proch;
  }
  else suiv=suiv->next;
 }
 suiv= NULL;
proch= NULL; parcour=NULL; newsui=NULL;
 free (newsui);free (proch);free(parcour);
}


/* vraisemblance Cox partielle avec tableau*/
void coxtablo()
{double val1,val2,val3,varest,det,pibeta,expos;
 int hh1,hh2,h1,h2,v1,v2,nused=0,i,j,nls;


 for (i=0;i<nall;i++) dmat2[i]=0;
 for (i=0;i<nall;i++) mdvs2[i]=0;
 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}


 nls=0;
 suiv = base;
 do
 {if (suiv->tnbhapo>0)
  {if (nls==0)
   {for (i=0;i<n;i++) tablo[nls][i]=0;
    tabpi[nls]=0;
   }
   else
   {for (i=0;i<n;i++) tablo[nls][i]=tablo[nls-1][i];
    tabpi[nls]=tabpi[nls-1];
   }
   val1=0;for (j=0;j<ajust;j++) val1+=effest[nbhest+j]*suiv->z[j];
   if (haplozero==0)
   {h1=suiv->hapest[0];h2=suiv->hapest[1];hh1=coding(h1);hh2=coding(h2);
    if (hh1>0) {val1+=effest[hh1];}
    if (hh2>0) {val1+=effest[hh2];}
    if (nbadd>0)
    {for (j=0;j<nbadd;j++)
     {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
      {val1+=effest[nbhest+ajust+j];}
     }
    }
    for (j=0;j<intercov;j++) {val1+=suiv->z[tabint[j][1]-1]*effest[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));}
   }
   expos=exp(val1);
   for (i=0;i<ajust;i++) {tablo[nls][nitptp[nbhest+i]]+=expos*suiv->z[i];}
   if (haplozero==0)
   {if ((hh1>0) && (itptp[hh1]==1)) {tablo[nls][nitptp[hh1]]+=expos;}
    if ((hh2>0) && (itptp[hh2]==1)) {tablo[nls][nitptp[hh2]]+=expos;}
    if (hypoth>0)
    {for (i=0;i<hypoth;i++)
     {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<=tabhypo[i][0]);
      v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<=tabhypo[i][0]);
      if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {tablo[nls][nitptp[hh1]]+=expos;}
      if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {tablo[nls][nitptp[hh2]]+=expos;}
     }
    }
    if (interor==1)
    {for (j=0;j<nbhypor;j++)
     {for (i=1;i<nbor[j];i++)
      {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
       {if (hh1==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;
	                               tablo[nls][nitptp[tabinter[j][i][0]-1]]+=expos;
	                	      }
	if (hh2==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;
	                               tablo[nls][nitptp[tabinter[j][i][0]-1]]+=expos;
	 			      }
       }
       else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
       {if (hh1==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;
                                       tablo[nls][nitptp[tabinter[j][0][0]-1]]-=expos;
                                      }
        if (hh2==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;
                                       tablo[nls][nitptp[tabinter[j][0][0]-1]]-=expos;
	 	                      }
       }
       else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
       {if (hh1==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;}
	if (hh2==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;}
       }
       else
       {if (hh1==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;
                                       tablo[nls][nitptp[tabinter[j][i][0]-1]]+=expos;
                                       tablo[nls][nitptp[tabinter[j][0][0]-1]]-=expos;
	       			      }
	if (hh2==tabinter[j][i][1]-1) {tablo[nls][nitptp[tabinter[j][0][1]-1]]+=expos;
                                       tablo[nls][nitptp[tabinter[j][i][0]-1]]+=expos;
  	                               tablo[nls][nitptp[tabinter[j][0][0]-1]]-=expos;
				      }
       }
      }
     }
    }
    for (i=0;i<nbadd;i++)
    {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
     {tablo[nls][nitptp[nbhest+ajust+i]]+=expos;}
    }
    for (i=0;i<intercov;i++)
    {tablo[nls][nitptp[nbhest+ajust+nbadd+i]]+=expos*suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
    }
   }
   tabpi[nls]+=expos;
   nls++;
   if (suiv->phen[0]==1.0)
   {for (i=0;i<n;i++) {dmat2[i]=tablo[nls-1][i];}
    for (i=0;i<n;i++) {dmat2[i]/=(-tabpi[nls-1]);}
    for (i=0;i<ajust;i++) {dmat2[nitptp[nbhest+i]]+=suiv->z[i];}
    if (haplozero==0)
    {hh1=coding(suiv->hapest[0]);
     hh2=coding(suiv->hapest[1]);
     if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=1;}
     if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=1;}
     if (hypoth>0)
     {for (i=0;i<hypoth;i++)
      {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<=tabhypo[i][0]);
       v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<=tabhypo[i][0]);
       if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=1;}
       if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=1;}
      }
     }
     if (interor==1)
     {for (j=0;j<nbhypor;j++)
      {for (i=1;i<nbor[j];i++)
       {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
				       }
         if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
     				       }
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	                               }
         if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	                               }
        }
        else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
         if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
        }
        else
	{if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	                               }
         if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
				       }
        }
       }
      }
     }
     for (i=0;i<nbadd;i++)
     {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
      {dmat2[nitptp[nbhest+ajust+i]]+=1;}
     }
     for (i=0;i<intercov;i++)
     {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }

    }
    /*ATTENTION L'INTERCEPT PEUT NE PAS ETRE CALCULE D'OU (N-1) */
    for (i=0;i<n-1;i++) {mdvs2[i]+=dmat2[i+1];}
    for (i=0;i<n-1;i++) {for (j=0;j<n-1;j++) {mdvd2[i][j]+=dmat2[i+1]*dmat2[j+1];}}

   }
  }
  suiv=suiv->next;
 } while ( (suiv!=NULL) && (suiv->next!=NULL)) ;
 suiv=NULL;




 /*ATTENTION L'INTERCEPT PEUT NE PAS ETRE CALCULE D'OU (N-1) */
 sysl(mdvd2,n-1);

 for (i=0;i<n-1;i++)
  {modif2[i]=0;
   for (j=0;j<n-1;j++) modif2[i]+=mdvd2[i][j]*mdvs2[j];
  }

 for (i=1;i<nall;i++) if (nitptp[i]>-1) {effest[i]+=modif2[nitptp[i]-1];}

 if (hypoth>0)
 {for (i=0;i<hypoth;i++)
  {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
   v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
   if (v1!=1) {effest[v2-1]=effest[v1-1];}
   else {effest[v2-1]=0;}
  }
 }
 if (interor==1)
 {for (j=0;j<nbhypor;j++)
  {for (i=1;i<nbor[j];i++)
   {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]-effest[tabinter[j][0][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1];
	}
	else
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1]-effest[tabinter[j][0][0]-1];
    }
   }
  }
 }
 if (hypint>0)
 {for (i=0;i<hypint;i++)
  {effest[nbhest+ajust+nbadd+tabhypint[i][1]-1]=effest[nbhest+ajust+nbadd+tabhypint[i][0]-1];
  }
 }

 for (i=0;i<n;i++) tempx+=(modif2[i]*modif2[i]);


}







/*Algorithme NR pour vraisemblance Cox partielle*/
void coxrubin(double **vintra)
{double val1,val2,val3,varest,det,pibeta,expos;
 int hh1,hh2,h1,h2,v1,v2,nused=0,i,j;

 dindividu *parcour;


 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}


 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {nused+=1-((suiv->tnbhapo>0) && (suiv->phen[0]==1.0));
  if ( (suiv->tnbhapo>0) && (suiv->phen[0]==1.0))
  {for (i=0;i<n;i++) {dmat2[i]=0;}
   parcour=base;
   pibeta=0;
   /* CALCUL DE Pi(beta)*/
   while ((parcour!=NULL) && (parcour->next!=NULL)  && (parcour->phen[1]>=suiv->phen[1]))
   {if ((parcour->tnbhapo>0))
    {val1=0;
     for (j=0;j<ajust;j++) val1+=effest[nbhest+j]*parcour->z[j];
     if (haplozero==0)
     {h1=parcour->hapest[0];      h2=parcour->hapest[1];
      hh1=coding(h1);             hh2=coding(h2);
      if (hh1>0) {val1+=effest[hh1];}
      if (hh2>0) {val1+=effest[hh2];}
      if (nbadd>0)
      {for (j=0;j<nbadd;j++)
       {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
        {val1+=effest[nbhest+ajust+j];}
       }
      }
      for (j=0;j<intercov;j++)
      {val1+=parcour->z[tabint[j][1]-1]*effest[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
      }
     }
     expos=exp(val1);
     for (i=0;i<ajust;i++) {dmat2[nitptp[nbhest+i]]+=expos*parcour->z[i];}
     if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=expos;}
      if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=expos;}
      if (hypoth>0)
      {for (i=0;i<hypoth;i++)
       {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<=tabhypo[i][0]);
        v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<=tabhypo[i][0]);
        if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=expos;}
	    if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=expos;}
	   }
      }
      if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (i=1;i<nbor[j];i++)
        {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	     {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
	                                     dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
	       							    }
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
	                                     dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
	 								    }
    	 }
	     else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
         {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
                                        }
          if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
	 	                                }
    	 }
	     else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
         {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;}
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;}
	     }
         else
	     {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
		      						    }
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
  	                                     dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
		 							    }
         }
        }
       }
      }
      for (i=0;i<nbadd;i++)
      {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
       {dmat2[nitptp[nbhest+ajust+i]]+=expos;}
      }
      for (i=0;i<intercov;i++)
      {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=expos*parcour->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
      }
	 }
     pibeta+=expos;
    }
	parcour=parcour->next;
   }/* FIN DE PARCOUR*/
    parcour=NULL;

   for (i=0;i<n;i++) {dmat2[i]/=-pibeta;}

   for (i=0;i<ajust;i++) {dmat2[nitptp[nbhest+i]]+=suiv->z[i];}
   if (haplozero==0)
   {hh1=coding(suiv->hapest[0]);
    hh2=coding(suiv->hapest[1]);
    if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=1;}
    if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=1;}
    if (hypoth>0)
    {for (i=0;i<hypoth;i++)
     {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<=tabhypo[i][0]);
      v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<=tabhypo[i][0]);
      if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=1;}
	  if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=1;}
	 }
    }
    if (interor==1)
    {for (j=0;j<nbhypor;j++)
     {for (i=1;i<nbor[j];i++)
      {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	   {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
	                                   dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
	     							  }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
     								  }
	    }
	   else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
       {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	 	                              }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	 	                              }
       }
       else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
       {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
       }
       else
	   {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	                                  }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
     	 							  }
       }
      }
     }
    }
    for (i=0;i<nbadd;i++)
    {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
     {dmat2[nitptp[nbhest+ajust+i]]+=1;}
    }
    for (i=0;i<intercov;i++)
    {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
    }
   }

  /*ATTENTION L'INTERCEPT PEUT NE PAS ETRE CALCULE D'OU (N-1) */
   for (i=0;i<n-1;i++) {for (j=0;j<n-1;j++) {mdvd2[i][j]+=dmat2[i+1]*dmat2[j+1];}}
  }
  suiv=suiv->next;
 }


 sysl(mdvd2,n-1);

 for (i=1;i<nall;i++)
	for (j=1;j<nall;j++)
	  if ((nitptp[i]>-1) && (nitptp[j]>-1))  vintra[i][j]+=mdvd2[nitptp[i]-1][nitptp[j]-1];

/* NEW DAVID*/
 free(parcour);parcour=NULL;

}


void breslow1(double *frqsem,double *lgoddsem,double *vtabres)
{double val1,v1,vrais,h1,h2,like,p1,p2,vraisfa=0,vv2,pib2=0,pib0=0,pib1=0,expos1,expos2;
 int hh1,hh2/*,idx*/,i,j;
 dindividu *parcour;


 for (i=0;i<3;i++) {tabres[i]=0;vtabres[i]=0;}


  suiv = base;
   while ((suiv!=NULL) && (suiv->next!=NULL))      /*ATTENTION AU PB EVENTUEL DES DONNES MANQUANTES*/
  {if ( (suiv->tnbhapo>0) && (suiv->phen[0]==1.0))
   {parcour=base;
    pib2=0;pib0=0;pib1=0;
    while ((parcour!=NULL) && (parcour->next!=NULL) && (parcour->phen[1]>=suiv->phen[1]))
    {if ( /*(parcour->phen[1]>=suiv->phen[1]) &&*/ (parcour->tnbhapo>0) )
	 {expos1=0;
	  pib0+=1;
	  val1=0;
	  for (j=0;j<ajust;j++) val1+=lgoddsem[nbhest+j]*parcour->z[j];
	  expos1=exp(val1);expos2=0.0;
	  if (haplozero==0)
	  {p2=0.0;
	   for (i=0;i<parcour->tnbhapo;i++)
       {h1=frqsem[parcour->idnb[i][0]];h2=frqsem[parcour->idnb[i][1]];
	    v1=0;p1=0.0;
        if ( (h1>0) && (h2>0))
        {p1=h1*h2*(2-(parcour->idnb[i][0]==parcour->idnb[i][1]));
	     hh1=coding(parcour->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
	     hh2=coding(parcour->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
         if (nbadd>0)
         {for (j=0;j<nbadd;j++)
          {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
           {v1+=lgoddsem[nbhest+ajust+j];}
          }
         }
         for (j=0;j<intercov;j++)
         {v1+=parcour->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
         }
	    }
		expos2+=p1*exp(v1);
	    p2+=p1;
	   }
	   expos2/=p2;
	  }
	  else {expos2=1.0;}
	  pib2+=expos1*expos2;
      pib1+=expos1;
	 }
     parcour=parcour->next;
	}
	p2=0.0;
	v1=0;for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
    expos2=0.0;
    tabres[1]=exp(v1);
	if (haplozero==0)
	{p2=0.0;
	 for (i=0;i<suiv->tnbhapo;i++)
     {h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
	  v1=0;p1=0.0;
      if ( (h1>0) && (h2>0))
      {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
	   hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
	   hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
       if (nbadd>0)
       {for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {v1+=lgoddsem[nbhest+ajust+j];}
        }
       }
       for (j=0;j<intercov;j++)
       {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
       }
	  }
	  expos2+=p1*exp(v1);
	  p2+=p1;
	 }
	 expos2/=p2;
	}
	else expos2=1.0;
	tabres[2]=tabres[1]*=expos2;

   /* BRESLOW SANS AUCUN EFFET*/
	tabres[0]=1/pib0;
	vtabres[0]-=log(tabres[0]);
	if (tabres[1]>0) {vtabres[1]-=log(tabres[1]/pib1);}
	if (tabres[2]>0) {vtabres[2]-=log(tabres[2]/pib2);}
	//if (vrais>0) {vraisfa -= log(vrais/pib2);}
   }
   suiv= suiv->next;
  }

  //tabres=NULL;free (tabres);
  //return(vtabres);
  parcour=NULL;
  free (parcour);


}

void likematchpair(double *frqsem,double *lgoddsem,double *vtabres)
{double val1,v1,vrais,h1,h2,like,p1,p2,vraisfa=0,vv2,pib2=0,pib0=0,pib1=0,expos1,expos2;
 int hh1,hh2/*,idx*/,i,j;
 dindividu *parcour;


 for (i=0;i<3;i++) {tabres[i]=0;vtabres[i]=0;}

  suiv = base;
   while ((suiv!=NULL) && (suiv->next!=NULL))      /*ATTENTION AU PB EVENTUEL DES DONNES MANQUANTES*/
  {if ( (suiv->tnbhapo>0) && (suiv->phen[0]==1.0))
   {parcour=base;
    pib2=0;pib0=0;pib1=0;
	while ((parcour!=NULL) && (parcour->next!=NULL) && (parcour->phen[1]<=suiv->phen[1]))
    {if ((parcour->tnbhapo>0)  && (parcour->phen[1]==suiv->phen[1]))
	 {expos1=0;
	  pib0+=1;
	  val1=0;
	  for (j=0;j<ajust;j++) val1+=lgoddsem[nbhest+j]*parcour->z[j];
	  expos1=exp(val1);expos2=0.0;
	  if (haplozero==0)
	  {p2=0.0;
	   for (i=0;i<parcour->tnbhapo;i++)
       {h1=frqsem[parcour->idnb[i][0]];h2=frqsem[parcour->idnb[i][1]];
	    v1=0;p1=0.0;
        if ( (h1>0) && (h2>0))
        {p1=h1*h2*(2-(parcour->idnb[i][0]==parcour->idnb[i][1]));
	     hh1=coding(parcour->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
	     hh2=coding(parcour->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
         if (nbadd>0)
         {for (j=0;j<nbadd;j++)
          {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
           {v1+=lgoddsem[nbhest+ajust+j];}
          }
         }
         for (j=0;j<intercov;j++)
         {v1+=parcour->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
         }
	    }
		expos2+=p1*exp(v1);
	    p2+=p1;
	   }
	   expos2/=p2;
	  }
	  else {expos2=1.0;}
	  pib2+=expos1*expos2;
      pib1+=expos1;
	 }
     parcour=parcour->next;
	}
	p2=0.0;
	v1=0;for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
    expos2=0.0;
    tabres[1]=exp(v1);
	if (haplozero==0)
	{p2=0.0;
	 for (i=0;i<suiv->tnbhapo;i++)
     {h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
	  v1=0;p1=0.0;
      if ( (h1>0) && (h2>0))
      {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
	   hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
	   hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
       if (nbadd>0)
       {for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {v1+=lgoddsem[nbhest+ajust+j];}
        }
       }
       for (j=0;j<intercov;j++)
       {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
       }
	  }
	  expos2+=p1*exp(v1);
	  p2+=p1;
	 }
	 expos2/=p2;
	}
	else expos2=1.0;
	tabres[2]=tabres[1]*=expos2;

   /* BRESLOW SANS AUCUN EFFET*/
	tabres[0]=1/pib0;
	vtabres[0]-=log(tabres[0]);
	if (tabres[1]>0) {vtabres[1]-=log(tabres[1]/pib1);}
	if (tabres[2]>0) {vtabres[2]-=log(tabres[2]/pib2);}
   }
   suiv= suiv->next;
  }
  //tabres=NULL;  free (tabres);
  parcour=NULL;
  free (parcour);


}



/*Algorithme NR pour vraisemblance Cox partielle*/
void matchpair()
{double val1,val2,val3,varest,det,pibeta,expos;
 int hh1,hh2,h1,h2,v1,v2,nused=0,i,j;


 dindividu *parcour;



 for (i=0;i<nall;i++) mdvs2[i]=0;
 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}


 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {nused+=1-((suiv->tnbhapo>0) && (suiv->phen[0]==1.0));
  if ( (suiv->tnbhapo>0) && (suiv->phen[0]==1.0))
  {for (i=0;i<n;i++) {dmat2[i]=0;}
   parcour=base;
   pibeta=0;
   /* CALCUL DE Pi(beta)*/
   while ((parcour!=NULL) && (parcour->next!=NULL)  && (parcour->phen[1]<=suiv->phen[1]))
   {if ((parcour->tnbhapo>0)  && (parcour->phen[1]==suiv->phen[1]))
	{
 	val1=0;
	 //val1=2*effest[0];
	 for (j=0;j<ajust;j++) val1+=effest[nbhest+j]*parcour->z[j];
	 if (haplozero==0)
     {h1=parcour->hapest[0];      h2=parcour->hapest[1];
	  hh1=coding(h1);             hh2=coding(h2);
      if (hh1>0) {val1+=effest[hh1];}
      if (hh2>0) {val1+=effest[hh2];}
      if (nbadd>0)
      {for (j=0;j<nbadd;j++)
       {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
        {val1+=effest[nbhest+ajust+j];}
       }
      }
      for (j=0;j<intercov;j++)
      {val1+=parcour->z[tabint[j][1]-1]*effest[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
      }
	 }
	 expos=exp(val1);



	 for (i=0;i<ajust;i++) {dmat2[nitptp[nbhest+i]]+=expos*parcour->z[i];}
     if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=expos;}
      if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=expos;}
      if (hypoth>0)
      {for (i=0;i<hypoth;i++)
       {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<=tabhypo[i][0]);
        v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<=tabhypo[i][0]);
        if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=expos;}
	    if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=expos;}
	   }
      }
      if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (i=1;i<nbor[j];i++)
        {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	     {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
	                                     dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
	       							    }
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
	                                     dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
	 								    }
    	 }
	     else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
         {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
                                        }
          if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
	 	                                }
    	 }
	     else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
         {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;}
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;}
	     }
         else
	     {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
		      						    }
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=expos;
                                         dmat2[nitptp[tabinter[j][i][0]-1]]+=expos;
  	                                     dmat2[nitptp[tabinter[j][0][0]-1]]-=expos;
		 							    }
         }
        }
       }
      }
      for (i=0;i<nbadd;i++)
      {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
       {dmat2[nitptp[nbhest+ajust+i]]+=expos;}
      }
      for (i=0;i<intercov;i++)
      {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=expos*parcour->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
      }
	 }
     pibeta+=expos;
    }
	parcour=parcour->next;
   }/* FIN DE PARCOUR*/


   for (i=0;i<n;i++) {dmat2[i]/=-pibeta;}



   for (i=0;i<ajust;i++) {dmat2[nitptp[nbhest+i]]+=suiv->z[i];}
   if (haplozero==0)
   {hh1=coding(suiv->hapest[0]);
    hh2=coding(suiv->hapest[1]);
    if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=1;}
    if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=1;}
    if (hypoth>0)
    {for (i=0;i<hypoth;i++)
     {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<=tabhypo[i][0]);
      v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<=tabhypo[i][0]);
      if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=1;}
	  if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=1;}
	 }
    }
    if (interor==1)
    {for (j=0;j<nbhypor;j++)
     {for (i=1;i<nbor[j];i++)
      {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	   {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
	                                   dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
	     							  }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
     								  }
	    }
	   else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
       {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	 	                              }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	 	                              }
       }
       else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
       {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
       }
       else
	   {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	                                  }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
     	 							  }
       }
      }
     }
    }
    for (i=0;i<nbadd;i++)
    {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
     {dmat2[nitptp[nbhest+ajust+i]]+=1;}
    }
    for (i=0;i<intercov;i++)
    {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
    }
   }

  /*ATTENTION L'INTERCEPT PEUT NE PAS ETRE CALCULE D'OU (N-1) */
   for (i=0;i<n-1;i++) {mdvs2[i]+=dmat2[i+1];}
   for (i=0;i<n-1;i++) {for (j=0;j<n-1;j++) {mdvd2[i][j]+=dmat2[i+1]*dmat2[j+1];}}
  }
  suiv=suiv->next;
 }




 /*ATTENTION L'INTERCEPT PEUT NE PAS ETRE CALCULE D'OU (N-1) */
 sysl(mdvd2,n-1);
 for (i=0;i<n-1;i++)
  {modif2[i]=0;
   for (j=0;j<n-1;j++) modif2[i]+=mdvd2[i][j]*mdvs2[j];
  }

 for (i=1;i<nall;i++) if (nitptp[i]>-1) {effest[i]+=modif2[nitptp[i]-1];}

 if (hypoth>0)
 {for (i=0;i<hypoth;i++)
  {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
   v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
   if (v1!=1) {effest[v2-1]=effest[v1-1];}
   else {effest[v2-1]=0;}
  }
 }
 if (interor==1)
 {for (j=0;j<nbhypor;j++)
  {for (i=1;i<nbor[j];i++)
   {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]-effest[tabinter[j][0][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1];
	}
	else
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1]-effest[tabinter[j][0][0]-1];
    }
   }
  }
 }
 if (hypint>0)
 {for (i=0;i<hypint;i++)
  {effest[nbhest+ajust+nbadd+tabhypint[i][1]-1]=effest[nbhest+ajust+nbadd+tabhypint[i][0]-1];
  }
 }

 for (i=0;i<n;i++) tempx+=(modif2[i]*modif2[i]);
 parcour=NULL;
  free (parcour);

}


void fishpair(double *frqsem, double *lgoddsem,matrixp sfisher)
{int i,j,k,hh1,hh2,vv2,vv1;
 double val1,val2,val3,h1,h2,like,likeli,v1,v2,p1,pibet,expos1,expos2,p2;
 double deriv[maxn],derivodd[maxn],derivfrq[maxn],tempod[maxn];
 int place2[maxn];
 double mfish[maxn][maxn];
 double help[maxn][maxn];

 double **tempfq;
 dindividu *parcour;

 printf("Running Variance Estimation\n");


 nnt=nbhest+n-1;


 for (i=0;i<maxn;i++) place2[i]=-1;
 j=0;for (i=0;i<nbhhypo;i++) if (frqsem[i]>0) {place2[i]=j;j++;}


 tempfq=(double **) malloc((size_t) (nbhest*sizeof(double *)));
 for (i=0;i<nbhest;i++) tempfq[i]=(double *) malloc((size_t) (2*sizeof(double)));


 //mfish=(double **) malloc((size_t) (nnt*sizeof(double *)));
 //for (i=0;i<nnt;i++) mfish[i]=(double *) malloc((size_t) (nnt*sizeof(double)));

 for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) mfish[i][j]=0;
 for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) sfisher[i][j]=0;

 for (i=0;i<nnt;i++) {deriv[i]=0;}
 for (i=0;i<n;i++) {derivodd[i]=0;tempod[i]=0;}
 for (i=0;i<nbhest;i++) {derivfrq[i]=0;tempfq[i][0]=0;tempfq[i][1]=0;}

 suiv=base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
 {if ((suiv->tnbhapo>0) && (suiv->phen[0]==1.0))
  {for (i=0;i<nbhest;i++) derivfrq[i]=0;
   for (i=0;i<n;i++) derivodd[i]=0;
   for (i=0;i<nnt;i++) deriv[i]=0;
   parcour = base;
   pibet=0;
   while ((parcour!=NULL) && (parcour->next!=NULL) && (parcour->phen[1]<=suiv->phen[1]))
   {if ((parcour->tnbhapo>0)  && (parcour->phen[1]==suiv->phen[1]))
     {
	 val1=0;
	 for (j=0;j<ajust;j++) val1+=lgoddsem[nbhest+j]*parcour->z[j];
     expos1=exp(val1);
     for (j=0;j<n;j++) {tempod[j]=0;}

	 for (j=0;j<nbhest;j++) {tempfq[j][0]=0;tempfq[j][1]=0;}
	 if (haplozero==1)
	 {for (j=0;j<ajust;j++) {derivodd[nitp[nbhest+j]]+=expos1*parcour->z[j];}
      pibet+=expos1;
	 }
	 else
     {p2=0.0;expos2=0.0;
	  for (i=0;i<parcour->tnbhapo;i++)
      {h1=frqsem[parcour->idnb[i][0]];h2=frqsem[parcour->idnb[i][1]];
	   p1=0.0;v1=0;
       if ( (h1>0) && (h2>0))
	   {p1=h1*h2*(2-(parcour->idnb[i][0]==parcour->idnb[i][1]));
        hh1=coding(parcour->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
        hh2=coding(parcour->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}


	    if (nbadd>0)
        {for (j=0;j<nbadd;j++)
         {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
          {v1+=lgoddsem[nbhest+ajust+j];}
         }
        }
		for (j=0;j<intercov;j++)
        {v1+=parcour->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
        }

        val1=exp(v1);
        if (parcour->idnb[i][0]==parcour->idnb[i][1])
	    {tempfq[place2[parcour->idnb[i][0]]][0]+=2*frqsem[parcour->idnb[i][0]]*val1;
         tempfq[place2[parcour->idnb[i][0]]][1]+=2*frqsem[parcour->idnb[i][0]];
        }
        else
	    {tempfq[place2[parcour->idnb[i][0]]][0]+=2*frqsem[parcour->idnb[i][1]]*val1;
         tempfq[place2[parcour->idnb[i][1]]][0]+=2*frqsem[parcour->idnb[i][0]]*val1;
	     tempfq[place2[parcour->idnb[i][0]]][1]+=2*frqsem[parcour->idnb[i][1]];
         tempfq[place2[parcour->idnb[i][1]]][1]+=2*frqsem[parcour->idnb[i][0]];
	    }
	    val1=p1*exp(v1);
        p2+=p1;
        expos2+=val1;

	    if ((hh1>0) && (itp[hh1]==1))  {tempod[nitp[hh1]]+=val1;}
	    if ((hh2>0) && (itp[hh2]==1))  {tempod[nitp[hh2]]+=val1;}
        if (hypoth>0)
		{for (j=0;j<hypoth;j++)
	     {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
		  v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
		  if ( (hh1==v2-1) && (nitp[hh1]>-1) ) {tempod[nitp[hh1]]+=val1;}
		  if ( (hh2==v2-1) && (nitp[hh2]>-1) ) {tempod[nitp[hh2]]+=val1;}
		 }
	    }
        if (interor==1)
        {for (j=0;j<nbhypor;j++)
	     {for (k=1;k<nbor[j];k++)
		  {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
		   {if (hh1==tabinter[j][k][1]-1)
		    {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][k][0]-1]]+=val1;
            }
		    if (hh2==tabinter[j][k][1]-1)
		    {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][k][0]-1]]+=val1;
		    }
		   }
		   else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
		   {if (hh1==tabinter[j][k][1]-1)
		    {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][0][0]-1]]-=val1;
            }
		    if (hh2==tabinter[j][k][1]-1)
		    {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][0][0]-1]]-=val1;
		    }
		   }
		   else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
		   {if (hh1==tabinter[j][k][1]-1) {tempod[nitp[tabinter[j][0][1]-1]]+=val1;}
		    if (hh2==tabinter[j][k][1]-1) {tempod[nitp[tabinter[j][0][1]-1]]+=val1;}
		   }
		   else
		   {if (hh1==tabinter[j][k][1]-1)
		    {tempod[nitp[tabinter[j][0][1]-1]]+=val1;tempod[nitp[tabinter[j][k][0]-1]]+=val1;
			 tempod[nitp[tabinter[j][0][0]-1]]-=val1;
		    }
		    if (hh2==tabinter[j][k][1]-1)
		    {tempod[nitp[tabinter[j][0][1]-1]]+=val1;tempod[nitp[tabinter[j][k][0]-1]]+=val1;
			 tempod[nitp[tabinter[j][0][0]-1]]-=val1;
		    }
		   }
		  }
         }
		}
        for (j=0;j<nbadd;j++)
		{if ( ((tadd[j][0]-1==hh1) &&  (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) &&  (tadd[j][1]-1==hh1)))
	     {tempod[nitp[nbhest+ajust+j]]+=val1;}
		}
		for (j=0;j<intercov;j++)
		{tempod[nitp[nbhest+ajust+nbadd+j]]+=val1*parcour->z[tabint[j][1]-1]*((hh1==(tabint[j][0]-1))+(hh2==(tabint[j][0]-1)));
	    }
       }
	  }
	  for (i=0;i<nbhest;i++) {derivfrq[i]+=(expos1*(tempfq[i][0]*p2-tempfq[i][1]*expos2)/(p2*p2));}
      //for (i=0;i<nbhest;i++) {derivfrq[i]+=expos1*(tempfq[i][0]/p2);}
	  for (i=0;i<n;i++)
      {if ((i<=nitp[nbhest-1]) || (i>nitp[nbhest+ajust-1])) {derivodd[i]+=tempod[i]*expos1/p2;}
      }
      expos2/=p2;
	  for (j=0;j<ajust;j++) {derivodd[nitp[nbhest+j]]+=expos1*expos2*parcour->z[j];}
	  pibet+=(expos1*expos2);
     }
	}
	parcour=parcour->next;
   }
   for (i=0;i<n;i++) derivodd[i]/=(-pibet);
   for (i=0;i<ajust;i++) derivodd[nitp[nbhest+i]]+=suiv->z[i];
   for (i=0;i<nbhest;i++) derivfrq[i]/=(-pibet);


   for (i=0;i<n;i++) {tempod[i]=0;}
   for (j=0;j<nbhest;j++) {tempfq[j][0]=0;tempfq[j][1]=0;}


   if (haplozero==0)
   {expos2=0.0;p2=0.0;
    for (i=0;i<suiv->tnbhapo;i++)
    {hh1=coding(suiv->idnb[i][0]);hh2=coding(suiv->idnb[i][1]);
     h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
     if ( (h1>0) && (h2>0))
	 {p1=0.0;v1=0;
      p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
      p2+=p1;
      if (hh1>0) {v1+=lgoddsem[hh1];}
	  if (hh2>0) {v1+=lgoddsem[hh2];}
      for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {v1+=lgoddsem[nbhest+ajust+j];}
      }
      for (j=0;j<intercov;j++)
      {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
      }
	  val1=exp(v1);
      if (suiv->idnb[i][0]==suiv->idnb[i][1])
	  {tempfq[place2[suiv->idnb[i][0]]][0]+=2*frqsem[suiv->idnb[i][0]]*val1;
       tempfq[place2[suiv->idnb[i][0]]][1]+=2*frqsem[suiv->idnb[i][0]];
	  }
      else
	  {tempfq[place2[suiv->idnb[i][0]]][0]+=2*frqsem[suiv->idnb[i][1]]*val1;
       tempfq[place2[suiv->idnb[i][1]]][0]+=2*frqsem[suiv->idnb[i][0]]*val1;
	   tempfq[place2[suiv->idnb[i][0]]][1]+=2*frqsem[suiv->idnb[i][1]];
       tempfq[place2[suiv->idnb[i][1]]][1]+=2*frqsem[suiv->idnb[i][0]];
	  }
	  val1=p1*exp(v1);
	  expos2+=val1;
      if ((hh1>0) && (itp[hh1]==1))  {tempod[nitp[hh1]]+=val1;}
	  if ((hh2>0) && (itp[hh2]==1))  {tempod[nitp[hh2]]+=val1;}
      if (hypoth>0)
	  {for (j=0;j<hypoth;j++)
	   {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
	    v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
	    if ( (hh1==v2-1) && (nitp[hh1]>-1) ) {tempod[nitp[hh1]]+=val1;}
	    if ( (hh2==v2-1) && (nitp[hh2]>-1) ) {tempod[nitp[hh2]]+=val1;}
	   }
	  }
      if (interor==1)
      {for (j=0;j<nbhypor;j++)
	   {for (k=1;k<nbor[j];k++)
	    {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
	     {if (hh1==tabinter[j][k][1]-1)
	      {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][k][0]-1]]+=val1;
          }
		  if (hh2==tabinter[j][k][1]-1)
		  {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][k][0]-1]]+=val1;
		  }
		 }
		 else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
		 {if (hh1==tabinter[j][k][1]-1)
		  {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][0][0]-1]]-=val1;
          }
		  if (hh2==tabinter[j][k][1]-1)
		  {tempod[nitp[tabinter[j][0][1]-1]]+=val1;    tempod[nitp[tabinter[j][0][0]-1]]-=val1;
		  }
		 }
		 else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
		 {if (hh1==tabinter[j][k][1]-1) {tempod[nitp[tabinter[j][0][1]-1]]+=val1;}
		  if (hh2==tabinter[j][k][1]-1) {tempod[nitp[tabinter[j][0][1]-1]]+=val1;}
		 }
		 else
		 {if (hh1==tabinter[j][k][1]-1)
		  {tempod[nitp[tabinter[j][0][1]-1]]+=val1;tempod[nitp[tabinter[j][k][0]-1]]+=val1;
		   tempod[nitp[tabinter[j][0][0]-1]]-=val1;
		  }
		  if (hh2==tabinter[j][k][1]-1)
		  {tempod[nitp[tabinter[j][0][1]-1]]+=val1;tempod[nitp[tabinter[j][k][0]-1]]+=val1;
		   tempod[nitp[tabinter[j][0][0]-1]]-=val1;
		  }
		 }
		}
       }
	  }
      for (j=0;j<nbadd;j++)
	  {if ( ((tadd[j][0]-1==hh1) &&  (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) &&  (tadd[j][1]-1==hh1)))
	   {tempod[nitp[nbhest+ajust+j]]+=val1;}
	  }
	  for (j=0;j<intercov;j++)
	  {tempod[nitp[nbhest+ajust+nbadd+j]]+=val1*suiv->z[tabint[j][1]-1]*((hh1==(tabint[j][0]-1))+(hh2==(tabint[j][0]-1)));
	  }
   	 }
	}
	for (i=0;i<n;i++)
    {if ((i<=nitp[nbhest-1]) || (i>nitp[nbhest+ajust-1])) {derivodd[i]+=(tempod[i]/expos2);}
    }
   for (i=0;i<nbhest;i++)  {derivfrq[i]+=((tempfq[i][0]/expos2)-(tempfq[i][1]/p2));}
   // for (i=0;i<nbhest;i++)  {derivfrq[i]+=(tempfq[i][0]/expos2);}
   }

   for (i=0;i<nbhest-1;i++) {deriv[i]=derivfrq[i+1];}
   for (i=0;i<n-1;i++) deriv[nbhest-1+i]=derivodd[i+1];
   for (i=0;i<(nnt-1);i++) for (j=0;j<(nnt-1);j++) mfish[i][j]+=(deriv[i]*deriv[j]);
  }
  suiv=suiv->next;
 }

  //help=(double **) malloc((size_t) ((n-1)*sizeof(double *)));
  //for (i=0;i<n-1;i++) help[i]=(double *) malloc((size_t) ((n-1)*sizeof(double)));

  for (i=0;i<n-1;i++) for (j=0;j<n-1;j++) help[i][j]=mfish[nbhest-1+i][nbhest-1+j];
  sysl(help,n-1);
  for (i=0;i<n-1;i++) for (j=0;j<n-1;j++) sfisher[nbhest-1+i][nbhest-1+j]=help[i][j];


  tempfq=NULL;  parcour=NULL;
  //mfish=NULL;free (mfish);
  free (parcour);free (tempfq);

  //help=NULL;free (help);
}


void pairedfish(double **vintra)
{int i,j,k,hh1,hh2,vv2,vv1;
 double val1,val2,val3,h1,h2,like,likeli,v1,v2,p1,pibet,expos1,expos2,p2;
 double tempod[maxn];
 dindividu *parcour;



 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}



 for (i=0;i<n;i++) {dmat2[i]=0;tempod[i]=0;}

 suiv=base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {if ((suiv->tnbhapo>0) && (suiv->phen[0]==1.0))
  {for (i=0;i<n;i++) dmat2[i]=0;
   parcour = base;
   pibet=0;
   while ((parcour!=NULL) && (parcour->next!=NULL) && (parcour->phen[1]<=suiv->phen[1]))
   {if ((parcour->tnbhapo>0)  && (parcour->phen[1]==suiv->phen[1]))
    {val1=0;
	 for (j=0;j<ajust;j++) val1+=effest[nbhest+j]*parcour->z[j];
     expos1=exp(val1);
     for (j=0;j<n;j++) {tempod[j]=0;}
     if (haplozero==1)
	 {for (j=0;j<ajust;j++) {dmat2[nitptp[nbhest+j]]+=expos1*parcour->z[j];}
      pibet+=expos1;
	 }
	 else
     {p2=0.0;expos2=0.0;v1=0.0;
	  hh1=coding(parcour->hapest[0]);if (hh1>0) {v1+=effest[hh1];}
      hh2=coding(parcour->hapest[1]);if (hh2>0) {v1+=effest[hh2];}
      for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {v1+=effest[nbhest+ajust+j];}
      }
      for (j=0;j<intercov;j++)
      {v1+=parcour->z[tabint[j][1]-1]*effest[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
      }
      val1=exp(v1); p2=1.0; expos2+=val1;
      if ((hh1>0) && (itptp[hh1]==1))  {tempod[nitptp[hh1]]+=val1;}
      if ((hh2>0) && (itptp[hh2]==1))  {tempod[nitptp[hh2]]+=val1;}

      for (j=0;j<hypoth;j++)
	  {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
	   v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
	   if ( (hh1==v2-1) && (nitptp[hh1]>-1) ) {tempod[nitptp[hh1]]+=val1;}
	   if ( (hh2==v2-1) && (nitptp[hh2]>-1) ) {tempod[nitptp[hh2]]+=val1;}
	  }
	  for (j=0;j<nbhypor;j++)
	  {for (k=1;k<nbor[j];k++)
	   {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
	    {if (hh1==tabinter[j][k][1]-1)
	     {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
         }
	     if (hh2==tabinter[j][k][1]-1)
	     {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
	     }
	    }
	    else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
	    {if (hh1==tabinter[j][k][1]-1)
	     {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
         }
	     if (hh2==tabinter[j][k][1]-1)
	     {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
	     }
	    }
	    else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
	    {if (hh1==tabinter[j][k][1]-1) {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;}
	     if (hh2==tabinter[j][k][1]-1) {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;}
	    }
	    else
	    {if (hh1==tabinter[j][k][1]-1)
	     {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
	      tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
	     }
	     if (hh2==tabinter[j][k][1]-1)
	     {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
	      tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
	     }
	    }
	   }
      }
	  for (j=0;j<nbadd;j++)
	  {if ( ((tadd[j][0]-1==hh1) &&  (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) &&  (tadd[j][1]-1==hh1)))
	   {tempod[nitptp[nbhest+ajust+j]]+=val1;}
	  }
	  for (j=0;j<intercov;j++)
	  {tempod[nitptp[nbhest+ajust+nbadd+j]]+=val1*parcour->z[tabint[j][1]-1]*((hh1==(tabint[j][0]-1))+(hh2==(tabint[j][0]-1)));
	  }

	  for (i=0;i<n;i++) /* verif contrainte*/
      {if ((i<=nitptp[nbhest-1]) || (i>nitptp[nbhest+ajust-1])) {dmat2[i]+=tempod[i]*expos1/p2;}
      }
      expos2/=p2;
	  for (j=0;j<ajust;j++) {dmat2[nitptp[nbhest+j]]+=expos1*expos2*parcour->z[j];}
	  pibet+=(expos1*expos2);
     }
	}
	parcour=parcour->next;
   }
   for (i=0;i<n;i++) dmat2[i]/=(-pibet);
   for (i=0;i<ajust;i++) dmat2[nitptp[nbhest+i]]+=suiv->z[i];

   for (i=0;i<n;i++) {tempod[i]=0;}

   if (haplozero==0)
   {expos2=0.0;
   	hh1=coding(suiv->hapest[0]);hh2=coding(suiv->hapest[1]);
    v1=0;p2=1;
	if (hh1>0) {v1+=effest[hh1];}
	if (hh2>0) {v1+=effest[hh2];}
    for (j=0;j<nbadd;j++)
    {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
     {v1+=effest[nbhest+ajust+j];}
    }
    for (j=0;j<intercov;j++)
    {v1+=suiv->z[tabint[j][1]-1]*effest[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
    }
	val1=exp(v1);
    expos2+=val1;
    if ((hh1>0) && (itptp[hh1]==1))  {tempod[nitptp[hh1]]+=val1;}
	if ((hh2>0) && (itptp[hh2]==1))  {tempod[nitptp[hh2]]+=val1;}
    for (j=0;j<hypoth;j++)
	{v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
	 v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
	 if ( (hh1==v2-1) && (nitptp[hh1]>-1) ) {tempod[nitptp[hh1]]+=val1;}
	 if ( (hh2==v2-1) && (nitptp[hh2]>-1) ) {tempod[nitptp[hh2]]+=val1;}
	}
	for (j=0;j<nbhypor;j++)
	{for (k=1;k<nbor[j];k++)
	 {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
	  {if (hh1==tabinter[j][k][1]-1)
	   {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;    tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
       }
	   if (hh2==tabinter[j][k][1]-1)
	   {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;    tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
	   }
	  }
	  else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
	  {if (hh1==tabinter[j][k][1]-1)
	   {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;    tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
       }
	   if (hh2==tabinter[j][k][1]-1)
	   {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;    tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
	   }
	  }
	  else if ((tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
	  {if (hh1==tabinter[j][k][1]-1) {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;}
	   if (hh2==tabinter[j][k][1]-1) {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;}
	  }
	  else
	  {if (hh1==tabinter[j][k][1]-1)
	   {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
	    tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
	   }
	   if (hh2==tabinter[j][k][1]-1)
	   {tempod[nitptp[tabinter[j][0][1]-1]]+=val1;tempod[nitptp[tabinter[j][k][0]-1]]+=val1;
	    tempod[nitptp[tabinter[j][0][0]-1]]-=val1;
	   }
	  }
	 }
    }
	for (j=0;j<nbadd;j++)
	{if ( ((tadd[j][0]-1==hh1) &&  (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) &&  (tadd[j][1]-1==hh1)))
	 {tempod[nitptp[nbhest+ajust+j]]+=val1;}
	}
	for (j=0;j<intercov;j++)
	{tempod[nitptp[nbhest+ajust+nbadd+j]]+=val1*suiv->z[tabint[j][1]-1]*((hh1==(tabint[j][0]-1))+(hh2==(tabint[j][0]-1)));
    }

    for (i=0;i<n;i++)
    {if ((i<=nitptp[nbhest-1]) || (i>nitptp[nbhest+ajust-1])) {dmat2[i]+=(tempod[i]/expos2);}
    }
  }
  for (i=0;i<n-1;i++) {for (j=0;j<n-1;j++) {mdvd2[i][j]+=dmat2[i+1]*dmat2[j+1];}}
  }
  suiv=suiv->next;
 }


 /*ATTENTION L'INTERCEPT PEUT NE PAS ETRE CALCULE D'OU (N-1) */
 sysl(mdvd2,n-1);

 for (i=1;i<nall;i++)
  for (j=1;j<nall;j++)
	  if ((nitptp[i]>-1) && (nitptp[j]>-1))  vintra[i][j]+=mdvd2[nitptp[i]-1][nitptp[j]-1];

 parcour=NULL;free (parcour);
}





void fisherscoring()  /* NEW ONE CORRIGEE POUR LE BUG DE PLP*/
{double val1,val2,val3,varest,det;
 int hh1,hh2,v1,v2,nused=0,i,j;
 char rep;


 for (i=0;i<n;i++) mdvs2[i]=0;
 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}



 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {if (suiv->tnbhapo>0)
  {for (i=0;i<n;i++) {dmat2[i]=0;}
   val1=suiv->phen[0];
   val2=2*effest[0];
   for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
   if ((chxt==1) && (offset==1)) val2+=suiv->z[ajust];
   if (haplozero==0)
   {hh1=coding(suiv->hapest[0]);if (hh1>0) {val2+=effest[hh1];}
    hh2=coding(suiv->hapest[1]);if (hh2>0) {val2+=effest[hh2];}
    for (i=0;i<nbadd;i++)
    {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
     {val2+=effest[nbhest+ajust+i];}
    }
    for (i=0;i<intercov;i++)
    {val2+=suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));}
   }
   if ((int) chxt==1)
   {val3=exp(val2)/(1+exp(val2));
    val1=-val3;varest=val1;
	dmat2[0]+=2*val1;
	for (i=0;i<ajust;i++)  {dmat2[nitptp[nbhest+i]]+=suiv->z[i]*val1;}

    if (haplozero==0)
    {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=val1;}
     if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=val1;}

     if (hypoth>0)
     {for (i=0;i<hypoth;i++)
       {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
        v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
        if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=val1;}
	    if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=val1;}
	   }
      }
      if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (i=1;i<nbor[j];i++)
        {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	     {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;
	                                   dmat2[nitptp[tabinter[j][i][0]-1]]+=val1;
	       							 }
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;
	                                     dmat2[nitptp[tabinter[j][i][0]-1]]+=val1;
	 								 }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=val1;
	 	                              }
         if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=val1;
	 	                              }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;}
	     if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;}
	    }
        else
	    {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=val1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=val1;
		      						 }
	     if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=val1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=val1;
  	                                    dmat2[nitptp[tabinter[j][0][0]-1]]-=val1;
		  							 }
        }
       }
      }
     }
     for (i=0;i<nbadd;i++)
     {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
      {dmat2[nitptp[nbhest+ajust+i]]+=val1;}
     }
     for (i=0;i<intercov;i++)
     {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=val1*suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }
    }
	if ((int) suiv->phen[0]==1)
	{dmat2[0]+=2;
     for (i=0;i<ajust;i++)  {dmat2[nitptp[nbhest+i]]+=suiv->z[i];}
     if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=1;}
      if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=1;}

      for (i=0;i<hypoth;i++)
      {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
       v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
       if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=1;}
	   if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=1;}
	  }
      if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (i=1;i<nbor[j];i++)
        {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	     {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
	                                   dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
	       							 }
	      if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
	                                     dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
	 								 }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	 	                              }
         if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
	 	                              }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
	     if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;}
	    }
        else
	    {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
		      						 }
	     if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=1;
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=1;
  	                                    dmat2[nitptp[tabinter[j][0][0]-1]]-=1;
		  							 }
        }
       }
      }
     }
     for (i=0;i<nbadd;i++)
     {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
      {dmat2[nitptp[nbhest+ajust+i]]+=1;}
     }
     for (i=0;i<intercov;i++)
     {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }
    }

	}

    for (i=0;i<n;i++) {mdvs2[i]+=dmat2[i];}
    for (i=0;i<n;i++) {for (j=0;j<n;j++) {mdvd2[i][j]+=dmat2[i]*dmat2[j];}}
  }
  else if ((int) chxt==2)
  {val1-=val2;varest=1;
   if (chxt==2) {dmat2[0]=2*varest;}
   for (i=0;i<ajust;i++)  dmat2[nitptp[nbhest+i]]+=dmat2[0]*suiv->z[i];
   if (haplozero==0)
   {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=dmat2[0];}
    if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=dmat2[0];}
    for (i=0;i<hypoth;i++)
     {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
      v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
      if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=dmat2[0];}
	  if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=dmat2[0];}
	 }
    for (j=0;j<nbhypor;j++)
    {for (i=1;i<nbor[j];i++)
     {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
      {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
	                                  dmat2[nitptp[tabinter[j][i][0]-1]]+=dmat2[0];
	       							 }
       if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
	                                   dmat2[nitptp[tabinter[j][i][0]-1]]+=dmat2[0];
	 								 }
      }
      else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
      {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
                                      dmat2[nitptp[tabinter[j][0][0]-1]]-=dmat2[0];
	 	                             }
       if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=dmat2[0];
	 	                             }
      }
      else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
      {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];}
       if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];}
      }
      else
	  {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=dmat2[0];
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=dmat2[0];
		      						 }
	   if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
                                       dmat2[nitptp[tabinter[j][i][0]-1]]+=dmat2[0];
  	                                   dmat2[nitptp[tabinter[j][0][0]-1]]-=dmat2[0];
		 							 }
      }
     }
    }
    for (i=0;i<nbadd;i++)
    {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
     {dmat2[nitptp[nbhest+ajust+i]]+=dmat2[0];}
    }
    for (i=0;i<intercov;i++) {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=dmat2[0]*suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));}
   }
   for (i=0;i<n;i++) {mdvs2[i]+=dmat2[i]*val1/varest;}
   for (i=0;i<n;i++) {for (j=0;j<n;j++) {mdvd2[i][j]+=dmat2[i]*dmat2[j]/varest;}}


  }
 }
  suiv=suiv->next;
 }


 sysl(mdvd2,n);
 for (i=0;i<n;i++)
  {modif2[i]=0;
   for (j=0;j<n;j++) modif2[i]+=mdvd2[i][j]*mdvs2[j];
  }

 j=0;for (i=0;i<nall;i++) if (nitptp[i]>-1) {effest[i]+=modif2[nitptp[i]];j++;}

 if (hypoth>0)
 {for (i=0;i<hypoth;i++)
  {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
   v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
   if (v1!=1) {effest[v2-1]=effest[v1-1];}
   else {effest[v2-1]=0;}
  }
 }
 if (interor==1)
 {for (j=0;j<nbhypor;j++)
  {for (i=1;i<nbor[j];i++)
   {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]-effest[tabinter[j][0][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1];
	}
	else
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1]-effest[tabinter[j][0][0]-1];
    }
   }
  }
 }
 if (hypint>0)
 {for (i=0;i<hypint;i++)
  {effest[nbhest+ajust+nbadd+tabhypint[i][1]-1]=effest[nbhest+ajust+nbadd+tabhypint[i][0]-1];
  }
 }

 if (chxt==2)
 {suiv = base;
  steres=0;
   while ((suiv!=NULL) && (suiv->next!=NULL))
  {if (suiv->tnbhapo>0)
   {val1=suiv->phen[0];val2=2*effest[0];
    for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
    if (haplozero==0)
    {hh1=coding(suiv->hapest[0]);if (hh1>0) {val2+=effest[hh1];}
     hh2=coding(suiv->hapest[1]);if (hh2>0) {val2+=effest[hh2];}

	 if (nbadd>0)
     {for (i=0;i<nbadd;i++)
      {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
       {val2+=effest[nbhest+ajust+i];}
      }
     }
     for (i=0;i<intercov;i++)
     {val2+=suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }
    }
    steres+=(val1-val2)*(val1-val2);
    nused+=1;
   }
   suiv=suiv->next;
  }
  ste=sqrt(steres/(nused-1));

 }

 for (i=0;i<n;i++) tempx+=(modif2[i]*modif2[i]);



}





void Xfisherscoring()
{double val1,val2,val3,varest,det;
 int hh1,hh2,v1,v2,nused=0,i,j;


 for (i=0;i<n;i++) mdvs2[i]=0;
 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}

 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {if (suiv->tnbhapo>0)
  {if ((int) suiv->z[0]==1)
   {for (i=0;i<n;i++) {dmat2[i]=0;}
    val1=suiv->phen[0];
    val2=effest[0];
    for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
    if (offset==1)  val2+=suiv->z[ajust];
	if (haplozero==0)
    {hh1=coding(suiv->hapest[0]);if (hh1>0) {val2+=0.5*effest[hh1];}
     hh2=coding(suiv->hapest[1]);if (hh2>0) {val2+=0.5*effest[hh2];}
 	 if (nbadd>0)
     {for (i=0;i<nbadd;i++)
      {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
       {val2+=effest[nbhest+ajust+i];}
      }
     }
     for (i=0;i<intercov;i++)
     {val2+=0.5*suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }
    }
    if (chxt==1) {val3=exp(val2)/(1+exp(val2));val1-=val3;varest=val3*(1-val3);}
    else if (chxt==2) {val1-=val2;varest=1;}
    dmat2[0]=varest;
    for (i=0;i<ajust;i++)  {dmat2[nitptp[nbhest+i]]+=dmat2[0]*suiv->z[i];}
    if (haplozero==0)
    {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=0.5*dmat2[0];}
     if ((hh2>0) && (itptp[hh2]==1)) {dmat2[nitptp[hh2]]+=0.5*dmat2[0];}
     if (hypoth>0)
     {for (i=0;i<hypoth;i++)
      {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
       v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
       if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=0.5*dmat2[0];}
	   if ( (hh2==v2-1) && (nitptp[hh2]>-1)) {dmat2[nitptp[hh2]]+=0.5*dmat2[0];}
	  }
     }
     if (interor==1)
     {for (j=0;j<nbhypor;j++)
      {for (i=1;i<nbor[j];i++)
       {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	    {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];
	                                    dmat2[nitptp[tabinter[j][i][0]-1]]+=0.5*dmat2[0];
	        						   }
	    if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];
	                                   dmat2[nitptp[tabinter[j][i][0]-1]]+=0.5*dmat2[0];
	 								  }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=0.5*dmat2[0];
	 	                             }
        if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=0.5*dmat2[0];
	 	                              }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];}
	     if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];}
	    }
        else
	    {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=0.5*dmat2[0];
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=0.5*dmat2[0];
		        					   }
	     if (hh2==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=0.5*dmat2[0];
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=0.5*dmat2[0];
  	                                    dmat2[nitptp[tabinter[j][0][0]-1]]-=0.5*dmat2[0];
		  							 }
        }
       }
      }
     }
     for (i=0;i<nbadd;i++)
     {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
      {dmat2[nitptp[nbhest+ajust+i]]+=dmat2[0];}
     }
     for (i=0;i<intercov;i++)
     {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=0.5*dmat2[0]*suiv->z[tabint[i][1]-1]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }
    }
   }
   else if ((int) suiv->z[0]==0)
   {for (i=0;i<n;i++) {dmat2[i]=0;}
    val1=suiv->phen[0];
    val2=effest[0];
    for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
	if (offset==1)  val2+=suiv->z[ajust];
	if (haplozero==0)
    {hh1=coding(suiv->hapest[0]);if (hh1>0) {val2+=effest[hh1];}
     for (i=0;i<intercov;i++)
     {val2+=suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*(hh1==tabint[i][0]-1);
     }
    }
    if (chxt==1) {val3=exp(val2)/(1+exp(val2));val1-=val3;varest=val3*(1-val3);}
    else if (chxt==2) {val1-=val2;varest=1;}
    dmat2[0]=varest;
    for (i=0;i<ajust;i++)  {dmat2[nitptp[nbhest+i]]+=dmat2[0]*suiv->z[i];}

	if (haplozero==0)
    {if ((hh1>0) && (itptp[hh1]==1)) {dmat2[nitptp[hh1]]+=dmat2[0];}
     if (hypoth>0)
     {for (i=0;i<hypoth;i++)
      {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
       v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
       if ( (hh1==v2-1) && (nitptp[hh1]>-1)) {dmat2[nitptp[hh1]]+=dmat2[0];}
	  }
     }
     if (interor==1)
     {for (j=0;j<nbhypor;j++)
      {for (i=1;i<nbor[j];i++)
       {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	    {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
	                                    dmat2[nitptp[tabinter[j][i][0]-1]]+=dmat2[0];
	        						   }
	    }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
                                       dmat2[nitptp[tabinter[j][0][0]-1]]-=dmat2[0];
	 	                             }
        }
	    else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
        {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];}
	    }
        else
	    {if (hh1==tabinter[j][i][1]-1) {dmat2[nitptp[tabinter[j][0][1]-1]]+=dmat2[0];
                                        dmat2[nitptp[tabinter[j][i][0]-1]]+=dmat2[0];
                                        dmat2[nitptp[tabinter[j][0][0]-1]]-=dmat2[0];
		        					   }
	    }
       }
      }
     }
     for (i=0;i<intercov;i++)
     {dmat2[nitptp[nbhest+ajust+nbadd+i]]+=dmat2[0]*suiv->z[tabint[i][1]-1]*(hh1==tabint[i][0]-1);
     }
    }
   }
   for (i=0;i<n;i++) {mdvs2[i]+=dmat2[i]*val1/varest;}
   for (i=0;i<n;i++) {for (j=0;j<n;j++) {mdvd2[i][j]+=dmat2[i]*dmat2[j]/varest;}}
  }
  suiv=suiv->next;
 }

 sysl(mdvd2,n);
 for (i=0;i<n;i++)
  {modif2[i]=0;
   for (j=0;j<n;j++) modif2[i]+=mdvd2[i][j]*mdvs2[j];
  }

 j=0;for (i=0;i<nall;i++) if (nitptp[i]>-1) {effest[i]+=modif2[nitptp[i]];j++;}

 if (hypoth>0)
 {for (i=0;i<hypoth;i++)
  {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
   v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
   if (v1!=1) {effest[v2-1]=effest[v1-1];}
   else {effest[v2-1]=0;}
  }
 }
 if (interor==1)
 {for (j=0;j<nbhypor;j++)
  {for (i=1;i<nbor[j];i++)
   {if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]-effest[tabinter[j][0][0]-1];
	}
	else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
    {effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1];
	}
	else
	{effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1]-effest[tabinter[j][0][0]-1];
    }
   }
  }
 }
 if (hypint>0)
 {for (i=0;i<hypint;i++)
  {effest[nbhest+ajust+nbadd+tabhypint[i][1]-1]=effest[nbhest+ajust+nbadd+tabhypint[i][0]-1];
  }
 }

 if (chxt==2)
 {suiv = base;
  steres=0;
  while ((suiv!=NULL) && (suiv->next!=NULL))
  {if (suiv->tnbhapo>0)
   {val1=suiv->phen[0];
    val2=effest[0];
	for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
    if (haplozero==0)
    {if ((int) suiv->z[0]==0)
     {hh1=coding(suiv->hapest[0]);if (hh1>0) {val2+=effest[hh1];}
      for (i=0;i<intercov;i++)
      {val2+=suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*(hh1==tabint[i][0]-1);
      }
	 }
	 else if ((int) suiv->z[0]==1)
	 {hh1=coding(suiv->hapest[0]);if (hh1>0) {val2+=0.5*effest[hh1];}
      hh2=coding(suiv->hapest[1]);if (hh2>0) {val2+=0.5*effest[hh2];}
      if (nbadd>0)
      {for (i=0;i<nbadd;i++)
       {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
        {val2+=effest[nbhest+ajust+i];}
       }
      }
      for (i=0;i<intercov;i++)
      {val2+=0.5*suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
      }
	 }
    }
    steres+=(val1-val2)*(val1-val2);
    nused+=1;
   }
   suiv=suiv->next;
  }
  ste=sqrt(steres/(nused-1));

 }

 for (i=0;i<n;i++) tempx+=(modif2[i]*modif2[i]);



}


void fishem(double *frqsem, double *lgoddsem, matrixp sfisher)
{int i,j,k,hh1,hh2,vv2,vv1;
 double val1,val2,val3,h1,h2,like,likeli,v1,v2,p1;
 double *deriv,*derivodd,*derivfrq;int *place;

 double fisher[maxn][maxn];

 printf("Running Variance Estimation\n");

 nnt=nbhest+n-1;




 place=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
 for (i=0;i<nbhhypo;i++) place[i]=-1;
 j=0;for (i=0;i<nbhhypo;i++) if (frqsem[i]>0) {place[i]=j;j++;}

 //fisher=(double **) malloc((size_t) (nnt*sizeof(double *)));
 //for (i=0;i<nnt;i++) fisher[i]=(double *) malloc((size_t) (nnt*sizeof(double)));
 for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) fisher[i][j]=0;

 deriv=(double *) malloc((size_t) (nnt*sizeof(double)));
 for (i=0;i<nnt;i++) deriv[i]=0;

 derivodd=(double *) malloc((size_t) (n*sizeof(double)));
 for (i=0;i<n;i++) derivodd[i]=0;

 derivfrq=(double *) malloc((size_t) (nbhest*sizeof(double)));
 for (i=0;i<nbhest;i++) derivfrq[i]=0;


 suiv=base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
 {for (i=0;i<nbhest;i++) derivfrq[i]=0;
  for (i=0;i<n;i++) derivodd[i]=0;
  for (i=0;i<nnt;i++) deriv[i]=0;

  /*CALCUL DE P(Y,G)*/
  val1=suiv->phen[0];
  likeli=0;
  for (i=0;i<suiv->tnbhapo;i++)
  {like=0;
   h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
   if ( (h1>0) && (h2>0))
   {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
    v1=2*lgoddsem[0];
    for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
    if ( (chxt==1)	&& (offset==1)) v1+=suiv->z[ajust];
	if (haplozero==0)
	{hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
     hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
     if (nbadd>0)
     {for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {v1+=lgoddsem[nbhest+ajust+j];}
      }
     }
     for (j=0;j<intercov;j++)
     {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
     }
    }
	if (chxt==1) {like=exp(v1*val1)/(1+exp(v1));}
    else if (chxt==2) {v2=-0.5*(val1-v1)*(val1-v1)/(ste*ste);like=exp(v2)/(ste*sqrt(2.0*pi));}
	if (suiv->idnb[i][0]==suiv->idnb[i][1])
    {derivfrq[place[suiv->idnb[i][0]]]+=2*like*frqsem[suiv->idnb[i][0]];}
    else
    {derivfrq[place[suiv->idnb[i][0]]]+=2*like*frqsem[suiv->idnb[i][1]];
     derivfrq[place[suiv->idnb[i][1]]]+=2*like*frqsem[suiv->idnb[i][0]];
    }
    if (chxt==1)
    {val2=like*(1-like)*p1*((val1==1)-(val1==0));}
    else if (chxt==2)  {val2=p1*like*(val1-v1)/(ste*ste);}
    derivodd[0]+=2*val2;
	for (k=0;k<ajust;k++) {derivodd[nitp[nbhest+k]]+=suiv->z[k]*val2;}
	if (haplozero==0)
	{if ((hh1>0) && (itp[hh1]==1)) {derivodd[nitp[hh1]]+=val2;}
	 if ((hh2>0) && (itp[hh2]==1)) {derivodd[nitp[hh2]]+=val2;}
	 if (hypoth>0)
	 {for (k=0;k<hypoth;k++)
      {vv1=tabhypo[k][0]*(tabhypo[k][0]<tabhypo[k][1])+tabhypo[k][1]*(tabhypo[k][1]<tabhypo[k][0]);
       vv2=tabhypo[k][1]*(tabhypo[k][0]<tabhypo[k][1])+tabhypo[k][0]*(tabhypo[k][1]<tabhypo[k][0]);
	   if ((hh1==vv2-1) && (nitp[hh1]>-1)) {derivodd[nitp[hh1]]+=val2;}
	   if ((hh2==vv2-1) && (nitp[hh2]>-1)) {derivodd[nitp[hh2]]+=val2;}
	  }
	 }
	 if (interor==1)
	 {for (j=0;j<nbhypor;j++)
      {for (k=1;k<nbor[j];k++)
       {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
	    {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
	                                    derivodd[nitp[tabinter[j][k][0]-1]]+=val2;
								 }
	     if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
	                                    derivodd[nitp[tabinter[j][k][0]-1]]+=val2;
	        						    }
  	    }
	    else if ( (tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
        {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
                                        derivodd[nitp[tabinter[j][0][0]-1]]-=val2;
		                                }
         if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
                                        derivodd[nitp[tabinter[j][0][0]-1]]-=val2;
		                                }
	    }
	    else if ( (tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
        {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;}
	     if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;}
	    }
        else
	    {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
                                        derivodd[nitp[tabinter[j][k][0]-1]]+=val2;
  		                                derivodd[nitp[tabinter[j][0][0]-1]]-=val2;
									    }
	     if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
                                        derivodd[nitp[tabinter[j][k][0]-1]]+=val2;
  		                                derivodd[nitp[tabinter[j][0][0]-1]]-=val2;
									    }
        }
       }
      }
	 }
     for (k=0;k<nbadd;k++)
     {if ( ((tadd[k][0]-1==hh1) && (tadd[k][1]-1==hh2)) || ((tadd[k][0]-1==hh2) && (tadd[k][1]-1==hh1)) )
      {derivodd[nitp[nbhest+ajust+k]]+=val2;}
     }
	 for (k=0;k<intercov;k++)
     {derivodd[nitp[nbhest+ajust+nbadd+k]]+=val2*suiv->z[tabint[k][1]-1]*((hh1==tabint[k][0]-1)+(hh2==tabint[k][0]-1));
     }
    }/* if haplozero*/
    likeli+=p1*like;
   }/* if */
  }/* for i*/
  /* ATTENTION A LA PREMIERE FREQUENCE D'OU LE -1 */
  if (likeli>0)
  {for (i=0;i<nbhest-1;i++) {deriv[i]=derivfrq[i+1]/likeli;}
   for (i=0;i<n;i++) {deriv[nbhest-1+i]=derivodd[i]/likeli;}
   for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) fisher[i][j]+=deriv[i]*deriv[j];
  }
  suiv=suiv->next;
 }/* while */
 printf("Inverting Variance Matrix....\n");
 sysl(fisher,nnt);
 for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) sfisher[i][j]=fisher[i][j];

 deriv=NULL;derivfrq=NULL;derivodd=NULL;place=NULL;
 free (deriv);free (derivfrq);free (derivodd);free(place);
 //fisher=NULL;free (fisher);
}

void Xfishem(double *frqsem, double *lgoddsem, matrixp sfisher)
{int i,j,k,hh1,hh2,vv2,vv1;
 double val1,val2,val3,h1,h2,like,likeli,v1,v2,p1;
 double *deriv,*derivodd,*derivfrq;int *place;

 double fisher[maxn][maxn];
 printf("Running Variance Estimation\n");

 nnt=nbhest+n-1;

 place=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
 for (i=0;i<nbhhypo;i++) place[i]=-1;
 j=0;for (i=0;i<nbhhypo;i++) if (frqsem[i]>0) {place[i]=j;j++;}

 //fisher=(double **) malloc((size_t) (nnt*sizeof(double *)));
 //for (i=0;i<nnt;i++) fisher[i]=(double *) malloc((size_t) (nnt*sizeof(double)));
 for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) fisher[i][j]=0;

 deriv=(double *) malloc((size_t) (nnt*sizeof(double)));
 for (i=0;i<nnt;i++) deriv[i]=0;

 derivodd=(double *) malloc((size_t) (n*sizeof(double)));
 for (i=0;i<n;i++) derivodd[i]=0;

 derivfrq=(double *) malloc((size_t) (nbhest*sizeof(double)));
 for (i=0;i<nbhest;i++) derivfrq[i]=0;


 suiv=base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
 {for (i=0;i<nbhest;i++) derivfrq[i]=0;
  for (i=0;i<n;i++) derivodd[i]=0;
  for (i=0;i<nnt;i++) deriv[i]=0;
  if ((int) suiv->z[0]==1)
  {/*CALCUL DE P(Y,G)*/
   val1=suiv->phen[0];
   likeli=0;
   for (i=0;i<suiv->tnbhapo;i++)
   {like=0;
    h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
    if ( (h1>0) && (h2>0))
    {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
     v1=lgoddsem[0];
     for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
 	 if ( (chxt==1)	&& (offset==1)) v1+=suiv->z[ajust];
	 if (haplozero==0)
	 {hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=0.5*lgoddsem[hh1];}
      hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=0.5*lgoddsem[hh2];}
      if (nbadd>0)
      {for (j=0;j<nbadd;j++)
       {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
        {v1+=lgoddsem[nbhest+ajust+j];}
       }
      }
      for (j=0;j<intercov;j++)
      {v1+=0.5*suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
      }
     }
	 if (chxt==1) {like=exp(v1*val1)/(1+exp(v1));}
     else if (chxt==2) {v2=-0.5*(val1-v1)*(val1-v1)/(ste*ste);like=exp(v2)/(ste*sqrt(2.0*pi));}
   	 if (suiv->idnb[i][0]==suiv->idnb[i][1])
     {derivfrq[place[suiv->idnb[i][0]]]+=2*like*frqsem[suiv->idnb[i][0]];}
     else
     {derivfrq[place[suiv->idnb[i][0]]]+=2*like*frqsem[suiv->idnb[i][1]];
      derivfrq[place[suiv->idnb[i][1]]]+=2*like*frqsem[suiv->idnb[i][0]];
     }
    /*Calcul des dérivees premieres par rapport aux parametres de regression */
     if (chxt==1)
     {val2=like*(1-like)*p1*((val1==1)-(val1==0));}
     else if (chxt==2)  {val2=p1*like*(val1-v1)/(ste*ste);}
     derivodd[0]+=val2;
	 for (k=0;k<ajust;k++) {derivodd[nitp[nbhest+k]]+=suiv->z[k]*val2;}
	 if (haplozero==0)
	 {if ((hh1>0) && (itp[hh1]==1)) {derivodd[nitp[hh1]]+=0.5*val2;}
	  if ((hh2>0) && (itp[hh2]==1)) {derivodd[nitp[hh2]]+=0.5*val2;}
	  if (hypoth>0)
	  {for (k=0;k<hypoth;k++)
       {vv1=tabhypo[k][0]*(tabhypo[k][0]<tabhypo[k][1])+tabhypo[k][1]*(tabhypo[k][1]<tabhypo[k][0]);
       vv2=tabhypo[k][1]*(tabhypo[k][0]<tabhypo[k][1])+tabhypo[k][0]*(tabhypo[k][1]<tabhypo[k][0]);
	   if ((hh1==vv2-1) && (nitp[hh1]>-1)) {derivodd[nitp[hh1]]+=0.5*val2;}
	   if ((hh2==vv2-1) && (nitp[hh2]>-1)) {derivodd[nitp[hh2]]+=0.5*val2;}
	   }
	  }
	  if (interor==1)
	  {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
	     {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;
	                                    derivodd[nitp[tabinter[j][k][0]-1]]+=0.5*val2;
								 }
	      if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;
	                                    derivodd[nitp[tabinter[j][k][0]-1]]+=0.5*val2;
	        						    }
  	     }
	     else if ( (tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
         {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;
                                        derivodd[nitp[tabinter[j][0][0]-1]]-=0.5*val2;
		                                }
          if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;
                                        derivodd[nitp[tabinter[j][0][0]-1]]-=0.5*val2;
		                                }
	     }
	     else if ( (tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
         {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;}
	      if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;}
	     }
         else
 	     {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;
                                        derivodd[nitp[tabinter[j][k][0]-1]]+=0.5*val2;
  		                                derivodd[nitp[tabinter[j][0][0]-1]]-=0.5*val2;
									    }
	      if (hh2==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=0.5*val2;
                                        derivodd[nitp[tabinter[j][k][0]-1]]+=0.5*val2;
  		                                derivodd[nitp[tabinter[j][0][0]-1]]-=0.5*val2;
									    }
         }
        }
       }
	  }
      for (k=0;k<nbadd;k++)
      {if ( ((tadd[k][0]-1==hh1) && (tadd[k][1]-1==hh2)) || ((tadd[k][0]-1==hh2) && (tadd[k][1]-1==hh1)) )
       {derivodd[nitp[nbhest+ajust+k]]+=val2;}
      }
	  for (k=0;k<intercov;k++)
      {derivodd[nitp[nbhest+ajust+nbadd+k]]+=0.5*val2*suiv->z[tabint[k][1]-1]*((hh1==tabint[k][0]-1)+(hh2==tabint[k][0]-1));
      }
     }/* if haplozero*/
     likeli+=p1*like;
    }/* if */
   }/* for i*/
  /* ATTENTION A LA PREMIERE FREQUENCE D'OU LE -1 */

  }
  else if ((int) suiv->z[0]==0)
  {
   val1=suiv->phen[0];
   likeli=0;
    for (i=0;i<suiv->tnbhapo;i++)
   {like=0;
    h1=frqsem[suiv->idnb[i][0]];
    if (h1>0)
    {p1=h1;
     v1=lgoddsem[0];
     for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
 	 if ( (chxt==1)	&& (offset==1)) v1+=suiv->z[ajust];
	 if (haplozero==0)
	 {hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
      for (j=0;j<intercov;j++)
      {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*(hh1==tabint[j][0]-1);
      }
     }
	 if (chxt==1) {like=exp(v1*val1)/(1+exp(v1));}
     else if (chxt==2) {v2=-0.5*(val1-v1)*(val1-v1)/(ste*ste);like=exp(v2)/(ste*sqrt(2.0*pi));}
   	 derivfrq[place[suiv->idnb[i][0]]]+=like;
     if (chxt==1) {val2=like*(1-like)*p1*((val1==1)-(val1==0));}
     else if (chxt==2)  {val2=p1*like*(val1-v1)/(ste*ste);}
     derivodd[0]+=val2;
	 for (k=0;k<ajust;k++) {derivodd[nitp[nbhest+k]]+=suiv->z[k]*val2;}
	 if (haplozero==0)
	 {if ((hh1>0) && (itp[hh1]==1)) {derivodd[nitp[hh1]]+=val2;}
	  if (hypoth>0)
	  {for (k=0;k<hypoth;k++)
       {vv1=tabhypo[k][0]*(tabhypo[k][0]<tabhypo[k][1])+tabhypo[k][1]*(tabhypo[k][1]<tabhypo[k][0]);
        vv2=tabhypo[k][1]*(tabhypo[k][0]<tabhypo[k][1])+tabhypo[k][0]*(tabhypo[k][1]<tabhypo[k][0]);
	    if ((hh1==vv2-1) && (nitp[hh1]>-1)) {derivodd[nitp[hh1]]+=val2;}
	   }
	  }
	  if (interor==1)
	  {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (tabinter[j][0][0]==1) && (tabinter[j][k][0]>1))
	     {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
	                                    derivodd[nitp[tabinter[j][k][0]-1]]+=val2;
								 }
	     }
	     else if ( (tabinter[j][k][0]==1) && (tabinter[j][0][0]>1))
         {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
                                        derivodd[nitp[tabinter[j][0][0]-1]]-=val2;
		                                }
         }
	     else if ( (tabinter[j][k][0]==1) && (tabinter[j][0][0]==1))
         {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;}
	     }
         else
 	     {if (hh1==tabinter[j][k][1]-1) {derivodd[nitp[tabinter[j][0][1]-1]]+=val2;
                                        derivodd[nitp[tabinter[j][k][0]-1]]+=val2;
  		                                derivodd[nitp[tabinter[j][0][0]-1]]-=val2;
									    }
	     }
        }
       }
	  }
      for (k=0;k<intercov;k++)
      {derivodd[nitp[nbhest+ajust+nbadd+k]]+=val2*suiv->z[tabint[k][1]-1]*(hh1==tabint[k][0]-1);
      }

     }/* if haplozero*/
     likeli+=p1*like;
    }/* if */
   }/* for i*/
  }
  if (likeli>0)
  {for (i=0;i<nbhest-1;i++) {deriv[i]=derivfrq[i+1]/likeli;}
   for (i=0;i<n;i++) {deriv[nbhest-1+i]=derivodd[i]/likeli;}
   for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) fisher[i][j]+=deriv[i]*deriv[j];
  }
  suiv=suiv->next;
 }/* while */
 printf("Inverting Variance Matrix....\n");
 sysl(fisher,nnt);
 for (i=0;i<nnt;i++) for (j=0;j<nnt;j++) sfisher[i][j]=fisher[i][j];



 deriv=NULL;derivfrq=NULL;derivodd=NULL;place=NULL;
 free (deriv);free (derivfrq);free (derivodd);free(place);
 //fisher=NULL;free (fisher);
}/* void*/

/*********************CALCUL DE LA VRAISEMBLANCE TOTALE ******************************************/
double likelihood(double *frqsem, double *lgoddsem)
{double val1,v1,vrais,h1,h2,like,p1,vraisfa=0,vv2;
 int hh1,hh2/*,idx*/,i,j;
  suiv = base;

  while ((suiv!=NULL) && (suiv->next!=NULL))       /*ATTENTION AU PB EVENTUEL DES DONNES MANQUANTES*/
  {val1=suiv->phen[0];
   if (suiv->tnbhapo>0)
   {vrais=0;
    for (i=0;i<suiv->tnbhapo;i++)
    {h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
     v1=2*lgoddsem[0];
     for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
     if ( (chxt==1)	&& (offset==1)) v1+=suiv->z[ajust];
	 if ( (h1>0) && (h2>0))
     {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
   	  if (haplozero==0)
	  {hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
       hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
       if (nbadd>0)
       {for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {v1+=lgoddsem[nbhest+ajust+j];}
        }
       }
       for (j=0;j<intercov;j++)
       {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
       }
      }
	  if (chxt==1)
 	  {like=exp(v1*val1)/(1+exp(v1));}
      else if (chxt==2)
      {vv2=-0.5*(val1-v1)*(val1-v1)/(ste*ste);
       like=exp(vv2)/(ste*sqrt(2.0*pi));
      }
	  vrais+=p1*like;
	 }
    }
    if (vrais>0) {vraisfa -= log(vrais);}
   }
   suiv= suiv->next;
  }
  return(-vraisfa);
}


/*********************CALCUL DE LA VRAISEMBLANCE TOTALE ******************************************/
double Xlikelihood(double *frqsem, double *lgoddsem)
{double val1,v1,vrais,h1,h2,like,p1,vraisfa=0,vv2;
 int hh1,hh2/*,idx*/,i,j;
  suiv = base;

  while ((suiv!=NULL) && (suiv->next!=NULL))       /*ATTENTION AU PB EVENTUEL DES DONNES MANQUANTES*/
  {val1=suiv->phen[0];
   vrais=0;
   if ((suiv->tnbhapo>0) && ((int) suiv->z[0]==1))
   {for (i=0;i<suiv->tnbhapo;i++)
    {h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
     v1=lgoddsem[0];
     for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
     if ( (chxt==1)	&& (offset==1)) v1+=suiv->z[ajust];
	 if ( (h1>0) && (h2>0))
     {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
   	  if (haplozero==0)
	  {hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=0.5*lgoddsem[hh1];}
       hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=0.5*lgoddsem[hh2];}
       if (nbadd>0)
       {for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {v1+=lgoddsem[nbhest+ajust+j];}
        }
       }
       for (j=0;j<intercov;j++)
       {v1+=0.5*suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
       }
      }
	  if (chxt==1)
 	  {like=exp(v1*val1)/(1+exp(v1));}
      else if (chxt==2)
      {vv2=-0.5*(val1-v1)*(val1-v1)/(ste*ste);
       like=exp(vv2)/(ste*sqrt(2.0*pi));
      }
	  vrais+=p1*like;
	 }
    }
   }
   else if ((suiv->tnbhapo==1) && ((int) suiv->z[0]==0))
   {h1=frqsem[suiv->idnb[0][0]];
    v1=lgoddsem[0];
    for (j=0;j<ajust;j++) v1+=lgoddsem[nbhest+j]*suiv->z[j];
    if ( (chxt==1)	&& (offset==1)) v1+=suiv->z[ajust];
	if (h1>0)
    {p1=h1;
   	 if (haplozero==0)
	 {hh1=coding(suiv->idnb[0][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
      for (j=0;j<intercov;j++)
      {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*(hh1==tabint[j][0]-1);
      }
     }
	 if (chxt==1) {like=exp(v1*val1)/(1+exp(v1));}
     else if (chxt==2)
	 {vv2=-0.5*(val1-v1)*(val1-v1)/(ste*ste);
      like=exp(vv2)/(ste*sqrt(2.0*pi));
     }
	 vrais+=p1*like;
	}
   }
   if (vrais>0) {vraisfa -= log(vrais);}
   suiv= suiv->next;
  }
  return(-vraisfa);
}






/*********************CALCUL DE LA VRAISEMBLANCE CONDITIONNELLE ******************************************/
double condlike(double *frqsem)
{double denume,pgeno,p1,h1,h2;
 int i;

  denume= 0.0;
  suiv = base;
   while ((suiv!=NULL) && (suiv->next!=NULL))
  {p1=0;pgeno=1.0;
   for (i=0;i<suiv->tnbhapo;i++)
   {h1=frqsem[suiv->idnb[i][0]];
    h2=frqsem[suiv->idnb[i][1]];
    if ( (h1>0) && (h2>0))
    {p1+=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
    }
   }
   pgeno*=p1;
   if (pgeno>0) {denume+=log(pgeno);}
   suiv = suiv->next;
  }
  return(denume);
}


double Xcondlike(double *frqsem)
{double denume,pgeno,p1,h1,h2;
 int i;

  denume= 0.0;
  suiv = base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
  {p1=0;pgeno=1.0;
   if ((int) suiv->z[0]==1)
   {for (i=0;i<suiv->tnbhapo;i++)
    {h1=frqsem[suiv->idnb[i][0]];
     h2=frqsem[suiv->idnb[i][1]];
     if ( (h1>0) && (h2>0))
     {p1+=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
     }
    }
   }
   else if ((int) suiv->z[0]==0)
   {if (suiv->tnbhapo>1) {printf("Male individuals must be unambiguous!\n");exit(0);}
    for (i=0;i<suiv->tnbhapo;i++)
    {h1=frqsem[suiv->idnb[i][0]];
     if (h1>0) p1+=h1;
    }
   }
   pgeno*=p1;
   if (pgeno>0) {denume+=log(pgeno);}
   suiv = suiv->next;
  }
  return(denume);
}


void phenomean(FILE *outf1,FILE *outf2,matrixp matse)
{dhaplotype *vect1;
 int i,j,k,kk;
 char l;
 double mean,valef,valf;


 fprintf(outf1,"\n\nExpected Phenotypic Mean [95%% CI] According to Estimated Haplotypes\n\n");
 fprintf(outf2,"<br><br>");
 fprintf(outf2,"<table align=center border=0  width=80%%>\n");
 fprintf(outf2,"<tr><td width=20%%> </td><td width=30%%> </td><td width=50%%> </td></tr>\n");
 fprintf(outf2,"<tr><td align=center colspan=""3"">Expected Phenotypic Mean [95%% CI] According to Estimated Haplotypes</td></tr>\n");
 fprintf(outf2,"<tr><td> </td></tr><tr><td> </td><td> </td></tr>\n");
 for (i=0;i<nbhest;i++)
 {vect1=tnbhbase;
  while ( (vect1!=NULL) && (fcoda2[vect1->numnew]!=numhap[i]))
  {vect1=vect1->down;}

  fprintf(outf2,"<tr><td align=center> ");
   if ((itp[i]==1) || ( (itp[i]==0) && (effest[i]!=0)) || ( (itp[i]==0) && (nitp[i]==-2)) )
  {for (kk=0;kk<nbloci;kk++)
   {l=letter[kk][0]*(vect1->listall[kk]==1)+letter[kk][1]*(vect1->listall[kk]==2);
    fprintf(outf1,"%c",l);fprintf(outf2,"%c",l);
   }
   mean=effest[0]+effest[i]*(i>0);
   fprintf(outf1,"\t%.5f ",mean);
   fprintf(outf2,"</td><td align=center>%.5f</td>",mean);
   if (i==0) {valef=sqrt(matse[nbhest-1][nbhest-1]);}
   else
   {valf=matse[nbhest-1][nbhest-1]+matse[nbhest-1+nitp[i]][nbhest-1+nitp[i]];
    valf+=2*matse[nbhest-1][nbhest-1+nitp[i]];
	valef=sqrt(valf);
   }
   fprintf(outf1,"[%.5f - %.5f]\n",mean-1.96*valef,mean+1.96*valef);
   fprintf(outf2,"<td align=left>[%.5f - %.5f]</td></tr>",mean-1.96*valef,mean+1.96*valef);
   vect1=NULL;
  }
 }

 fprintf(outf2,"</table>\n");
 /* NEW DAVID*/
 free (vect1);vect1=NULL;
 


}

double residuel(double *frqsem, double *lgoddsem)
{double residu=0,val1,v1,h1,h2,p1;
 int hh1,hh2,i,j,nuse=0;
  suiv = base;
   while ((suiv!=NULL) && (suiv->next!=NULL))
  {val1=suiv->phen[0]-2*lgoddsem[0];
   for (j=0;j<ajust;j++) val1-=lgoddsem[nbhest+j]*suiv->z[j];
   if (suiv->tnbhapo>0)
   {nuse+=1;
    for (i=0;i<suiv->tnbhapo;i++)
    {h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
     p1=0;
 	 if ( (h1>0) && (h2>0))
     {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
      v1=0;
      if (haplozero==0)
	  {hh1=coding(suiv->idnb[i][0]);if (hh1>0) {v1+=lgoddsem[hh1];}
       hh2=coding(suiv->idnb[i][1]);if (hh2>0) {v1+=lgoddsem[hh2];}
       if (nbadd>0)
       {for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {v1+=lgoddsem[nbhest+ajust+j];}
        }
       }
       for (j=0;j<intercov;j++)
       {v1+=suiv->z[tabint[j][1]-1]*lgoddsem[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
       }
      }
	  val1-=p1*v1;
     }
    }
   	residu+=val1*val1;
   }
   suiv= suiv->next;
  }
  residu/=(nuse-1);

  return(residu);
}



double somdelai()
{double somd=0.0;
 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {if (suiv->tnbhapo>0)
  // && (suiv->phen[0]==1.0))
  {somd+=suiv->phen[1];
  }
  suiv=suiv->next;
 }
 return(somd);
}








//fin fonc12.c



int lecture( char *fileName,
	     int pMaxvarfic,
	     int pNbloci,
	     int *pIdloci,
	     int pLdmatrix,
	     int pMsdata,
	     int pR2,
	     int pChxt,
	     int pNum0,
	     int pIdtime,
	     int pOffset,
	     int pIdoffset,
	     int pAjust,
	     int *pNumajust,
	     int pXlnk,
	     int pNumsx);



/**********************************Lecture du fichier de donnees*****************************/
int lecture( char *fileName, 	// 1 ok
	     int pMaxvarfic,	// 2 ok
	     int pNbloci,	// 3 ok
	     int *pIdloci,	// 4 ok
	     int pLdmatrix,	// 5 ok
	     int pMsdata,	// 6 ok
	     int pR2,		// 7 ok
	     int pChxt,		// 8 ok
	     int pNum0,		// 9 ok
	     int pIdtime,	// 10 ok
	     int pOffset,	// 16
	     int pIdoffset,	// 11 ok
	     int pAjust,	// 12 ok
	     int *pNumajust,	// 13 ok
	     int pXlnk,		// 14
	     int pNumsx)	// 15
{
	char flec[30],rep,cha;
 	int maxvarfic,mxvfic=0,i,idp,idx,numh,numsx,num0,dtmq,a=0,cpt_tab=0,l;
 	int *idloc,iref,idtime,idwgt;
 	double res;
 	vectgen geno;
 	FILE *fildo /* *fillet*/;
 	int *taille,*idlet;
 	char **tempo;

 	int kel = 0;
 	printf("No loci : %d\n",pNbloci);
 	for( kel = 0 ; kel < pNbloci ; kel++)
	{
	 	printf("loc[%d] = %d", kel+1, pIdloci[kel]);
 	}

 	for( kel = 0 ; kel < pAjust ; kel++)
 	{
		printf("pNumajust[%d] = %d", kel+1, pNumajust[kel]);
 	}

 	printf("DataFile name :\n");

//-------------------- 1 -----------------------//
 	strcpy(flec,fileName);
//-------------------- fin 1 -----------------------//

 	if ((fildo=fopen(flec,"r"))==NULL)
 	{
		return -1;
 	}

//-------------------- 2 -----------------------//
 	maxvarfic = pMaxvarfic;
//-------------------- fin 2 -----------------------//

 	printf("Number of di-allelic loci to be studied :\n");

//-------------------- 3 -----------------------//
 	nbloci = pNbloci;
//-------------------- fin 3 -----------------------//

 	idloc =(int *)malloc((size_t) (nbloci*sizeof(int)));

//-------------------- 4 -----------------------//
 	for (i=0;i<nbloci;i++)
 	{
	 	idloc[i] = pIdloci[i];
 	}
//-------------------- fin 4 -----------------------//

 	alfreq=(double *) malloc((size_t) (nbloci*sizeof(double)));
 	tnbhbase=(dhaplotype *)malloc((size_t) (sizeof(dhaplotype)));

 	for (i=1;i<=nbloci;i++)
 	{
 		idx=ipow(2,nbloci-i);
 		iref=1;
 		idp=1;
 		numh=0;
  		tnbhnew=tnbhbase;
  		do
  		{
  			tnbhnew->listall[i-1]=iref;
   			if (i==1)
   			{
	   			tnbhnew->numnew=numh;
	   			tnbhnew->present=0;
              			if (numh==ipow(2,nbloci)-1) {tnbhnew->down=NULL;}
				else {tnbhnew->down=(dhaplotype *) malloc((size_t) (sizeof(dhaplotype)));}
			}
   			idp++;
   			if (idp>idx)
   			{
   				if (iref==1) {iref=2;}
                 		else	     {iref=1;}
				idp=1;
			}
   			numh++;
   			tnbhnew=tnbhnew->down;
  		}
  		while (numh<ipow(2,nbloci));
  		tnbhnew=NULL;
 	}
//-------------------- 5 -----------------------//
 	ldmatrix = pLdmatrix;
//-------------------- fin 5 -----------------------//

//-------------------- 6 -----------------------//
 	msdata = pMsdata;
//-------------------- fin 6 -----------------------//

//-------------------- 7 -----------------------//
 	rsq = pR2;
//-------------------- fin 7 -----------------------//

//-------------------- 8 -----------------------//
 	chxt = pChxt;
//-------------------- fin 8 -----------------------//

//-------------------- 14 -----------------------//
	xlnk = pXlnk;
//-------------------- fin 14 -----------------------//

//-------------------- 15 -----------------------//
	numsx = pNumsx;
//-------------------- fin 15 -----------------------//

//-------------------- 16 -----------------------//
	offset = pOffset;
//-------------------- fin 16 -----------------------//

/* cas improbable !!
 	if ((chxt==0) && (xlnk==1))
 	{
 		ajust=1;
 		numajust[0]=numsx;
 	}
*/

 	if (chxt>0)
 	{
//-------------------- 9 -----------------------//
		num0 = pNum0;
//-------------------- fin 9 -----------------------//

//-------------------- 10 -----------------------//
  		if ((chxt==3) || (chxt==4))
  		{
 			idtime = pIdtime;
  		}
//-------------------- fin 10 -----------------------//

// GT 04-10-06 ce mecanisme est implementé dans GraficT..
		/*
		if (xlnk==0)
  		{
//-------------------- 12 -----------------------//
  			ajust=pAjust;
//-------------------- fin 12 -----------------------//
	 		if(ajust > 0)
	 		{
//-------------------- 13 -----------------------//
				for (i=0;i<ajust;i++)
				{
					numajust[i] = pNumajust[i];
				}
//-------------------- fin 13 -----------------------//
	 		}
  		}
  		else if (xlnk==1)
  		{

//-------------------- 12 -----------------------//
  			ajust = pAjust + 1;
//-------------------- fin 12 -----------------------//

  			numajust[0]=numsx;
  			if(ajust > 1)
	 		{
//-------------------- 13 -----------------------//
				for (i=0;i<ajust;i++)
				{
					numajust[i+1] = pNumajust[i];
				}
//-------------------- fin 13 -----------------------//
	 		}
  		}
		*/
//-------------------- 12 -----------------------//
  		ajust=pAjust;
//-------------------- fin 12 -----------------------//
	 	if(ajust > 0)
	 	{
//-------------------- 13 -----------------------//
			for (i=0;i<ajust;i++)
			{
				numajust[i] = pNumajust[i];
			}
//-------------------- fin 13 -----------------------//
		}
// Fin GT-1- 04-10-06
	}
	if (chxt==1)
  	{
  		offset = pOffset;
  		idoffset = 0;
  		if (offset == 1)
  		{
//-------------------- 11 -----------------------//
  			idoffset = pIdoffset;
//-------------------- fin 11 -----------------------//
		}
    		if (ajust+offset>maxcov)
    		{
    			printf("With such an offset and so many covaraites, the analysis can not be performed.\n");
     			printf("Please contact D.Tregouet for further information\n");
     			exit(0);
    		}
  	}

  	do
  	{
  		cha=fgetc(fildo);
   		if ( (cha=='\t') || (cha==32) || (cha==';'))
   		{
   			mxvfic++;
   		}
  	}while (cha!='\n');
  	mxvfic++;

 	rewind (fildo);
 	if (mxvfic!=maxvarfic)
 	{
 		printf("Error in number of variables in Datafile\n");
  		exit(0);
 	}

 	cpt_tab=0;
  	taille =(int *)malloc((size_t) (10 * maxvarfic*sizeof(int)));
  	tempo =(char **) malloc((size_t) (10 * maxvarfic*sizeof(char *)));
  	for (i=0;i<maxvarfic;i++) tempo[i] = (char *) malloc((size_t) (30*sizeof(char)));
  	idlet =(int *)malloc((size_t) (nbloci*sizeof(int)));
  	for (i=0;i<nbloci;i++) idlet[i]=0;

 	l=0;
 	base =(dindividu *)malloc((size_t) (sizeof(dindividu)));
 	suiv=base;
 	while (!feof(fildo))
 	{
 		/* on lit le caractere */
 		cha=fgetc(fildo);
		/* caractere normal */
  		if (cha!='*' && cha!='\n' && cha!='\t' && cha!='â' && cha!=';' && cha!=32)
  		{
  			if (a<30)
   			{
   				/* on le stocke */
   				flec[a]=cha;
   				/* on termine la chaine pour ne pas avoir d erreur */
    				flec[a+1]='\0';
   			}
   			/* on incremente la compteur de string pour le caractere suivant */
   			a++;
  		}
  		/* on rencontre une nouvelle variable */
  		else if ((cha=='\t') || (cha==';') || (cha==32))
  		{
  			/* on stocke tout le double recupere */
  			res=atof(flec);
  			/* tempo[cpt_tab]=res;*/
   			strncpy(tempo[cpt_tab],flec,a);
   			taille[cpt_tab]=a;
   			cpt_tab++;
   			a=0;
  		}
  		/* on rencontre un nouvel individu */
  		else if (cha=='\n')
  		{
  			/* on stocke tout le double recupere */
  			res= atof(flec);
  			/* tempo[cpt_tab]=res;*/
   			strncpy(tempo[cpt_tab],flec,a);
   			taille[cpt_tab]=a;
   			a=0;
   			flec[0]='\0';
   			cpt_tab=0;
   			dtmq=0;
   			for (i=0;i<nbloci;i++)
   			{/*flec[0]='\0';*/
    				strncpy(flec,tempo[idloc[i]-1],taille[idloc[i]-1]);
   				flec[taille[idloc[i]-1]]='\0';
				if (flec[0]=='0')
    				{
    					geno[i][0]=0;
    					geno[i][1]=0;
    					dtmq+=1;
    				}
    				else if (flec[0]==flec[1])
				{
					if (idlet[i]==2)
     					{
     						if (flec[0]==letter[i][0])
     						{
     							geno[i][0]=1;
     							geno[i][1]=1;
     						}
      						else
      						{
      							geno[i][0]=2;
      							geno[i][1]=2;
      						}
	 				}
     					else if (idlet[i]==1)
	 				{
	 					if (flec[0]==letter[i][0])
	 					{
	 						geno[i][0]=1;
	 						geno[i][1]=1;
	 					}
      						else
      						{
      							letter[i][1]=flec[0];
      							geno[i][0]=2;
      							geno[i][1]=2;
      							idlet[i]=2;
      						}
     					}
	 				else if (idlet[i]==0)
	 				{
	 					letter[i][0]=flec[0];
	 					geno[i][0]=1;
	 					geno[i][1]=1;
	 					idlet[i]=1;
	 				}
   				}
				else if (flec[0]!=flec[1])
				{
					geno[i][0]=1;
					geno[i][1]=2;
	 				if (idlet[i]==0)
	 				{
	 					letter[i][0]=flec[0];
	 					letter[i][1]=flec[1];
	 				}
	 				else if (idlet[i]==1)
		  			{
		  				if (letter[i][0]==flec[0])
		  				{
		  					letter[i][1]=flec[1];
		  				}
	       					else letter[i][1]=flec[0];
		  			}
	 				idlet[i]=2;
    				}

  			}
   			for (i=0;i<ajust;i++)
   			{
   				strncpy(flec,tempo[numajust[i]-1],taille[numajust[i]-1]);
				flec[taille[numajust[i]-1]]='\0';
				res=atof(flec);
    				suiv->z[i]=res;
   			}
   			if (chxt>0)
   			{
   				strncpy(flec,tempo[num0-1],taille[num0-1]);
				flec[taille[num0-1]]='\0';
				res=atof(flec);
    				suiv->phen[0]=res;
    				suiv->phen[1]=-1;
				if ( (chxt==1) && (offset==1))
				{
					strncpy(flec,tempo[idoffset-1],taille[idoffset-1]);
     					flec[taille[idoffset-1]]='\0';
	 				res=atof(flec);
	 				suiv->z[ajust]=res;
    				}
				if ( (chxt==3) || (chxt==4))
    				{
    					strncpy(flec,tempo[idtime-1],taille[idtime-1]);
     					flec[taille[idtime-1]]='\0';
	 				res=atof(flec);
	 				suiv->phen[1]=res;
    				}
   			}
   			for (i=0;i<nbloci;i++)
   			{
   				suiv->marq[i][0]=geno[i][0];
   				suiv->marq[i][1]=geno[i][1];
   			}
   			suiv->nblm=dtmq;
   			suiv->next=(dindividu *)malloc((size_t) (sizeof(dindividu)));	/* on alloue pour un individu de plus */
   			suiv=suiv->next;
                        suiv->next=NULL;					   		/* on decale la liste d un cran */
  		}
 	}
 	suiv=NULL;

 	if (fildo !=NULL) fclose(fildo); fildo=NULL;

  	idlet=NULL;
  	taille=NULL;
  	tempo=NULL;
  	idloc=NULL;
 	free (idloc);
 	free (tempo);
 	free (taille);
 	free (idlet);

 	return 0;
}

/*****************Determine les haplotypes de chaque individu: corps de la procedure***********/
void determhapo()
{int i ,nls=0;
 nbtotused=0;
 vectgen xg;
 printf("Running identification of haplotypes....\n");
 if (msdata==1)
 {suiv=base;
  printf("(It can take quite a long time since you are dealing with missing data...\n");
 while ((suiv!=NULL) && (suiv->next!=NULL))
  {for (i=0;i<nbloci;i++) {xg[i][0]=suiv->marq[i][0];xg[i][1]=suiv->marq[i][1];}
   suiv->tnbhapo=0;
   nbhapo1(xg);
   if (suiv->tnbhapo>0) nbtotused+=1;
   suiv=suiv->next;
  }
 }
 else
 { suiv=base;
  do
  {for (i=0;i<nbloci;i++) {xg[i][0]=suiv->marq[i][0];xg[i][1]=suiv->marq[i][1];}
     suiv->tnbhapo=0;
    nbhapo0(xg);
    if (suiv->tnbhapo>0) nbtotused+=1;
    suiv=suiv->next;
  }
  while ((suiv!=NULL) && (suiv->next!=NULL));

 }
suiv=NULL;
}


/******Determine les haplotypes de chaque individu: sous procedure avec donnees manquantes*****/
void nbhapo1(vectgen xg)
{int i,tnbh1;
 short keepon,h0,h00,l0,l00;
 tnbh1=0;
 if (suiv->nblm<nbloci-locmq+1)
 {vect1=tnbhbase;
  while (vect1!=NULL)
  {vect2=vect1;
   while (vect2!=NULL)
   {i=0;
    do
    {h0=xg[i][0];h00=xg[i][1];
     l0=vect1->listall[i];l00=vect2->listall[i];
     if ( (((h0==l0) || (h0==0)) && ((h00==l00) || (h00==0)))
	   || (((h0==l00) || (h0==0)) && ((h00==l0) || (h00==0))) )
     {i+=1;keepon=1;}
     else {keepon=0;}
    }
    while ((keepon==1) && (i!=nbloci));
    if ((keepon==1) && (i==nbloci))
    {tnbh1+=1;}
    vect2=vect2->down;
   }
   vect1=vect1->down;
  }
 }
 suiv->idnb=(int **)malloc((size_t) (tnbh1*sizeof(int *)));
 for (i=0;i<tnbh1;i++) suiv->idnb[i]=(int *)malloc((size_t) (2*sizeof(int)));
 if (tnbh1>maxhapair) {maxhapair=tnbh1;}
 suiv->tnbhapo=tnbh1;
 tnbh1=0;
 if (suiv->nblm<nbloci-locmq+1)
 {vect1=tnbhbase;
  while (vect1!=NULL)
  {vect2=vect1;
   while (vect2!=NULL)
   {i=0;
    do
    {h0=xg[i][0];h00=xg[i][1];
     l0=vect1->listall[i];l00=vect2->listall[i];
     if ( (((h0==l0) || (h0==0)) && ((h00==l00) || (h00==0)))
	  || (((h0==l00) || (h0==0)) && ((h00==l0) || (h00==0))) )
     {i+=1;keepon=1;}
     else {keepon=0;}
    }
    while ((keepon==1) && (i!=nbloci));
    if ((keepon==1) && (i==nbloci))
    {tnbh1+=1;
     suiv->idnb[tnbh1-1][0]=vect1->numnew;
     suiv->idnb[tnbh1-1][1]=vect2->numnew;
    }
    vect2=vect2->down;
   }
   vect1=vect1->down;
  }
 }
 for (i=0;i<tnbh1;i++)
   {fcoda1[suiv->idnb[i][0]]=1;fcoda1[suiv->idnb[i][1]]=1;}
}



/******Determine les haplotypes de chaque individu: sous procedure sans donnees manquantes*****/
void nbhapo0(vectgen xg)
{int i,nbhet,numh,h0,h00,backha,idp,idx,tnbh1;
 int hetero[maxloc];
 short iref;char rep;
 tnbh1=0;nbhet=0;
 if (suiv->nblm==0)
 {for (i=0;i<nbloci;i++)
  {if (xg[i][0]!=xg[i][1]) {nbhet+=1;}
  }
  if (nbhet<2) {tnbh1=1;}
  else {tnbh1=ipow(2,nbhet-1);}
 }
 if (tnbh1>0)
 {suiv->idnb=(int **)malloc((size_t) (tnbh1*sizeof(int *)));
  for (i=0;i<tnbh1;i++) suiv->idnb[i]=(int *)malloc((size_t) (2*sizeof(int)));
 }

 if (tnbh1>maxhapair) {maxhapair=tnbh1;}
 suiv->tnbhapo=tnbh1;




 nbhet=0;
 if (suiv->nblm==0)
 {for (i=0;i<nbloci;i++)
  {if (xg[i][0]!=xg[i][1]) {hetero[nbhet]=i+1;nbhet+=1;}

  }
  if ( (nbhet>0) && (xlnk==1) &&  ((int) suiv->z[0]==0))
  {

   printf("One male individual heterozygotes at %d loci\n",nbhet);
   printf("In this current version, the individual must be deleted\n");
   exit(0);
   printf("Do you want to continue (y/n) ?");
   scanf("%c",&rep);
   if ((rep=='N')  || (rep=='n')) {exit(0);}
  }
  if (nbhet<2)
  {h0=0;h00=0;
   for (i=1;i<=nbloci;i++)
   {h0+=(xg[i-1][0]==2)*ipow(2,nbloci-i);
    h00+=(xg[i-1][1]==2)*ipow(2,nbloci-i);
   }
   suiv->idnb[0][0]=h0;
   suiv->idnb[0][1]=h00;
   fcoda1[h0]=1;fcoda1[h00]=1;
  }
  else
  {backha=0;
   for (i=1;i<=nbloci;i++)
   {if (xg[i-1][0]==xg[i-1][1]) {backha+=(xg[i-1][0]==2)*ipow(2,nbloci-i);}
   }
   for (i=0;i<tnbh1;i++)
   {suiv->idnb[i][0]=backha;
  	suiv->idnb[i][1]=backha+ipow(2,nbloci-hetero[0]);
   }
   for (i=1;i<=nbhet-1;i++)
   {numh=0;idx=ipow(2,nbhet-1-i);iref=1;idp=1;
    do
	{suiv->idnb[numh][0]+=(iref==2)*ipow(2,nbloci-hetero[i]);
	 suiv->idnb[numh][1]+=(iref==1)*ipow(2,nbloci-hetero[i]);
	 idp++;
	 if (idp>idx) {if (iref==1) {iref=2;}
	               else         {iref=1;}
				   idp=1;
			      }
     numh++;
	}while (numh<tnbh1);
   }
   for (i=0;i<tnbh1;i++)
   {fcoda1[suiv->idnb[i][0]]=1;fcoda1[suiv->idnb[i][1]]=1;}
  }
 }
}




/**********************Calcul des frequences alleliques et haplotypiques sans LD***************/
void initfreq(double *hw)
{int i,l;
 int **ngeno,*nl,**ngenoh,**ngenof,nh,nf;
 double fi1,fi2;
 ngeno=(int **) malloc ((size_t) (nbloci*sizeof(int *)));
 for (i=0;i<nbloci;i++) ngeno[i]=(int *) malloc((size_t) (3*sizeof(int)));
 nl=(int *) malloc ((size_t) (nbloci*sizeof(int)));
 for (i=0;i<nbloci;i++) {nl[i]=0;for (l=0;l<3;l++) ngeno[i][l]=0;}



 for (i=0;i<nbloci;i++) alfreq[i]=0;
 suiv=base;
 l=0;
 if (xlnk==0)
 {if (msdata==0)
  { while ((suiv!=NULL) && (suiv->next!=NULL))
   {if (suiv->nblm==0)
    {l++;
     for (i=0;i<nbloci;i++)
	  {alfreq[i]+=(suiv->marq[i][0]==1)+(suiv->marq[i][1]==1);
       ngeno[i][0]+=( (suiv->marq[i][0]==1) && (suiv->marq[i][1]==1) );
       ngeno[i][2]+=( (suiv->marq[i][0]==2) && (suiv->marq[i][1]==2) );
       ngeno[i][1]+=( ((suiv->marq[i][0]==1) && (suiv->marq[i][1]==2)) || ((suiv->marq[i][0]==2) && (suiv->marq[i][1]==1)));
	  }
    }
    suiv=suiv->next;
   }
   for (i=0;i<nbloci;i++) alfreq[i]/=2*l;
   for (i=0;i<nbloci;i++) {hw[i]=0;}


   for (i=0;i<nbloci;i++)
   {fi1=alfreq[i];fi2=1-fi1;
    hw[i]+=pow((ngeno[i][0]-fi1*fi1*l),2)/(fi1*fi1*l);
    hw[i]+=pow((ngeno[i][1]-2*fi1*fi2*l),2)/(2*fi1*fi2*l);
    hw[i]+=pow((ngeno[i][2]-fi2*fi2*l),2)/(fi2*fi2*l);
   }
  }
  else if (msdata==1)
  {suiv=base;
   for (i=0;i<nbloci;i++) nl[i]=0;
   while ((suiv!=NULL) && (suiv->next!=NULL))
   {if (suiv->nblm<nbloci-locmq+1)
    {for (i=0;i<nbloci;i++)
     {if ( (suiv->marq[i][0]==1) || (suiv->marq[i][1]==2))
      {nl[i]+=1;
	   alfreq[i]+=(suiv->marq[i][0]==1)+(suiv->marq[i][1]==1);
	   ngeno[i][0]+=( (suiv->marq[i][0]==1) && (suiv->marq[i][1]==1) );
       ngeno[i][2]+=( (suiv->marq[i][0]==2) && (suiv->marq[i][1]==2) );
       ngeno[i][1]+=( ((suiv->marq[i][0]==1) && (suiv->marq[i][1]==2)) || ((suiv->marq[i][0]==2) && (suiv->marq[i][1]==1)));
	  }
	 }
    }
    suiv=suiv->next;
   }
   for (i=0;i<nbloci;i++) alfreq[i]/=(2*nl[i]);
   for (i=0;i<nbloci;i++) {hw[i]=0;}
   for (i=0;i<nbloci;i++)
   {fi1=alfreq[i];fi2=1-fi1;
    hw[i]+=pow((ngeno[i][0]-fi1*fi1*nl[i]),2)/(fi1*fi1*nl[i]);
    hw[i]+=pow((ngeno[i][1]-2*fi1*fi2*nl[i]),2)/(2*fi1*fi2*nl[i]);
    hw[i]+=pow((ngeno[i][2]-fi2*fi2*nl[i]),2)/(fi2*fi2*nl[i]);
   }
  }
 }
 else if (xlnk==1)
 {if (msdata==1)
  {printf("The combined missing data and Xlinked options cannot be handled \n");exit(0);
  }

  ngenoh=(int **) malloc ((size_t) (nbloci*sizeof(int *)));
  ngenof=(int **) malloc ((size_t) (nbloci*sizeof(int *)));
  for (i=0;i<nbloci;i++)
  {ngenoh[i]=(int *) malloc((size_t) (3*sizeof(int)));
   ngenof[i]=(int *) malloc((size_t) (3*sizeof(int)));
  }
  for (i=0;i<nbloci;i++) {for (l=0;l<3;l++) {ngenoh[i][l]=0;ngenof[i][l]=0;}}

  nf=0;nh=0;
  l=0;
  suiv=base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
  {if (suiv->nblm==0)
   {if ((int) suiv->z[0]==1)
    {l+=2;
	 nf++;
	 for (i=0;i<nbloci;i++)
	 {alfreq[i]+=(suiv->marq[i][0]==1)+(suiv->marq[i][1]==1);
      ngenof[i][0]+=( (suiv->marq[i][0]==1) && (suiv->marq[i][1]==1));
	  ngenof[i][2]+=( (suiv->marq[i][0]==2) && (suiv->marq[i][1]==2) );
	  ngenof[i][1]+=( ((suiv->marq[i][0]==1) && (suiv->marq[i][1]==2)) || ((suiv->marq[i][0]==2) && (suiv->marq[i][1]==1)));
	 }
    }
	else if ((int) suiv->z[0]==0)
    {l+=1;
     nh++;
	 for (i=0;i<nbloci;i++)
	 {alfreq[i]+=(suiv->marq[i][0]==1);
      ngenoh[i][0]+=(suiv->marq[i][0]==1) ;
	  ngenoh[i][2]+=(suiv->marq[i][0]==2) ;
	 }
    }
   }
   suiv=suiv->next;
  }

  for (i=0;i<nbloci;i++)
  {alfreq[i]/=l;
   hw[i]=0;
   fi1=alfreq[i];fi2=1-fi1;
   hw[i]+=pow((ngenof[i][0]-fi1*fi1*nbhf[1][0]),2)/(fi1*fi1*nbhf[1][0]);
   hw[i]+=pow((ngenoh[i][0]-fi1*nbhf[0][0]),2)/(fi1*nbhf[0][0]);
   hw[i]+=pow((ngenof[i][1]-2*fi1*fi2*nbhf[1][0]),2)/(2*fi1*fi2*nbhf[1][0]);
   hw[i]+=pow((ngenof[i][2]-fi2*fi2*nbhf[1][0]),2)/(fi2*fi2*nbhf[1][0]);
   hw[i]+=pow((ngenoh[i][2]-fi2*nbhf[0][0]),2)/(fi2*nbhf[0][0]);
  }
 }

 freqest=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
 tempfreq=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
 if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
 {freqdist=(double **) malloc((size_t) (nbhhypo*sizeof(double *)));
  tempdist=(double **) malloc((size_t) (nbhhypo*sizeof(double *)));
  for (i=0;i<nbhhypo;i++)
  {freqdist[i]=(double *) malloc((size_t) (2*sizeof(double)));
   tempdist[i]=(double *) malloc((size_t) (2*sizeof(double)));
  }
 }
 if (chxt==5)
 {freqdist=(double **) malloc((size_t) (nbhhypo*sizeof(double *)));
  tempdist=(double **) malloc((size_t) (nbhhypo*sizeof(double *)));
  for (i=0;i<nbhhypo;i++)
  {freqdist[i]=(double *) malloc((size_t) (nbcatego*sizeof(double)));
   tempdist[i]=(double *) malloc((size_t) (nbcatego*sizeof(double)));
  }
 }

 tnbhnew=tnbhbase;
 while (tnbhnew!=NULL)
 {if (tnbhnew->present==1)
  {tnbhnew->frqle=1.0;
   for (l=0;l<nbloci;l++)
   {tnbhnew->frqle*=alfreq[l]*(tnbhnew->listall[l]==1)+(1-alfreq[l])*(tnbhnew->listall[l]==2);
   }
   i=fcoda2[tnbhnew->numnew];
   freqest[i]=tnbhnew->frqle;
  }
  tnbhnew=tnbhnew->down;
 }
 ngeno=NULL;nl=NULL;ngenoh=NULL;ngenof=NULL;
// New DAVID
 free ((int **) ngeno);free ((int *) nl);free ((int **) ngenof);free ((int **) ngenoh);

}

/*************************Lecture des contraintes sur les effets *********************/
void lecteffe()
{
	char rep;int i,j,v1,v2;double vv;
 	short chgt;
 	FILE *parfile;
	FILE *readParam;
 	haplozero=1;

 	readParam = fopen("paramData.thi","r");

 	printf("Do you want to estimate haplotypic effect (y/n) ?\n");
 	do
 	{
 		fscanf(readParam,"%c",&rep);
 	}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

 	if ((rep=='y') || (rep=='Y'))
 	{
 		haplozero=0;
 	}
  	if ((haplozero==0) && (chxt>0))
  	{
  		hypoth=0;
  		printf("Do you want to test specific hypothesis on haplotypic effects (y/n) ?\n");
  		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

  		if ((rep=='y') || (rep=='Y'))
  		{
  			printf("How many constraints of equality do you wish to consider ?\n");
   			fscanf(readParam,"%d",&hypoth);
   			tabhypo=(int **)malloc((size_t) (hypoth*sizeof(int *)));
   			for (i=0;i<hypoth;i++)
   			{
   				tabhypo[i]=(int *)malloc((size_t) (2*sizeof(int)));
   			}
   			for (i=0;i<hypoth;i++)
   			{
   				printf("Hypothesis %d: enter the number of the two haplotypes whose effect must set to be equal\n",i+1);
    				fscanf(readParam,"%d %d",&v1,&v2);
				tabhypo[i][0]=v1;
				tabhypo[i][1]=v2;
   			}
  		}
  		interor=0;
  		printf("Do you want to test for the homogeneity of some allelic effects  (y/n) ?\n");
  		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

  		if ((rep=='y') || (rep=='Y'))
  		{
  			interor=1;
   			printf("How many different hypotheses do you wish to test ?\n");
   			do
   			{
   				fscanf(readParam,"%d",&nbhypor);
   			}while (nbhypor>maxhypor);

   			tabinter=(int ***)malloc((size_t) (nbhypor*sizeof(int **)));
   			nbor=(int *)malloc((size_t) (nbhypor*sizeof(int)));

   			for (j=0;j<nbhypor;j++)
   			{
   				chgt=0;
    				printf("Hypothesis %d: How many ORs or Differences to be tested ?\n",j+1);
				do
				{
					fscanf(readParam,"%d",&v1);
					nbor[j]=v1;
				}while (nbor[j]<2);

				tabinter[j]=(int **)malloc((size_t) (nbor[j]*sizeof(int *)));
				for (i=0;i<nbor[j];i++)
				{
					tabinter[j][i]=(int *)malloc((size_t) (2*sizeof(int)));
				}
				for (i=0;i<nbor[j];i++)
				{
					printf("Enter the number of the two haplotypes involved in OR/Difference %d\n",i+1);
     					fscanf(readParam,"%d %d",&v1,&v2);
	 				tabinter[j][i][0]=v1;
	 				tabinter[j][i][1]=v2;
	 				if (v2==1)
	 				{
	 					chgt=1;
	 				}
    				}
				if (chgt==1)
				{
					for (i=0;i<nbor[j];i++)
	 				{
	 					v1=tabinter[j][i][0];
	  					tabinter[j][i][0]=tabinter[j][i][1];
	  					tabinter[j][i][1]=v1;
     					}
				}
   			}
  		}
  		nbadd=0;
  		printf("Do you want to test the non-additivity of some haplotypic effects (y/n)\n");
  		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

  		if ((rep=='y') || (rep=='Y'))
  		{
  			printf("How many hypothesis to be tested ?\n");
   			do
   			{
   				fscanf(readParam,"%d",&nbadd);
   			}while (nbadd>maxadd);

   			tadd=(int **)malloc((size_t) (nbadd*sizeof(int *)));
   			for (j=0;j<nbadd;j++)
   			{
   				tadd[j]=(int *)malloc((size_t) (2*sizeof(int)));
				printf("Hypothesis %d: Enter the two haplotypic effects assumed to interact \n",j+1);
    				fscanf(readParam,"%d %d",&v1,&v2);
				tadd[j][0]=v1;
				tadd[j][1]=v2;
   			}
  		}
  		if (ajust>0)
 		{
 			printf("Do you want to test haplotype * environment interactions (y/n) ?");
   			do
   			{
   				fscanf(readParam,"%c",&rep);
   			}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));
   			if ((rep=='y') || (rep=='Y'))
   			{
   				printf("How many interactions to be tested ?\n");
    				do
				{
					fscanf(readParam,"%d",&intercov);
				}while (intercov>maxcov);

				tabint=(int **)malloc((size_t) (intercov*sizeof(int *)));
				for (j=0;j<intercov;j++)
				{
					tabint[j]=(int *)malloc((size_t) (2*sizeof(int)));
	 				printf("Interaction %d: Enter first the haplotype  and then the covariate number.\n",j+1);
     					fscanf(readParam,"%d %d",&v1,&v2);
	 				tabint[j][0]=v1;
	 				tabint[j][1]=v2;
    				}
    				hypint=0;
    				/*printf("Do you want to set equality constraints on some interactions (y/n) ?\n");
    				do
    				{
    					fscanf(readParam,"%c",&rep);
    				}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

    				if ((rep=='y') || (rep=='Y'))
    				{
    					printf("How many constraints ?\n");
     					do
     					{
     						fscanf(readParam,"%d",&hypint);
     					}while (hypint>maxhypint);

	 				tabhypint=(int **)malloc((size_t) (hypint*sizeof(int *)));

	 				for (j=0;j<hypint;j++)
	 				{
	 					tabhypint[j]=(int *)malloc((size_t) (2*sizeof(int)));
     						printf("hypothesis %d: Enter the two interaction numbers set to be equal.\n",j+1);
	  					fscanf(readParam,"%d %d",&v1,&v2);
      						tabhypint[j][0]=v1;
      						tabhypint[j][1]=v2;
	 				}
				}*/
   			}
  		}
 	}

	parfile=NULL;
  	if ((parfile=fopen("para.txt","r"))==NULL)
  	{
	  	printf("Error in Para file reading....\n");
	  	exit(0);
  	}
  	// Lecture du nombre d'entrées dan sle fichier para.txt
  	fscanf(parfile,"%d\n",&nbhest);

  	nall=nbhest+nbadd+intercov+ajust;
  	numhap=(int *)malloc((size_t) (nbhest*sizeof(int)));

  	for (i=0;i<nbhhypo;i++)
  	{
  		freqest[i]=0;
  	}
  	for (i=0;i<nbhest;i++)
  	{
  		// lecture de des entrées (4 parametres)
  		fscanf(parfile,"%d %lf %d %lf\n",&v2,&vv,&itp[i],&effest[i]);
		numhap[i]=fcoda2[v2];
   		freqest[numhap[i]]=vv;
	}

 	vv=0;
 	for (i=0;i<nbhhypo;i++)
 	{
 		vv+=freqest[i];
 	}

  	for (i=0;i<nbhhypo;i++)
  	{
  		freqest[i]/=vv;
  	}
  	vv=0;
  	for (i=0;i<nbhhypo;i++)
  	{
  		vv+=freqest[i];
  	}
  	for (i=nbhest;i<nall;i++)
  	{
  		itp[i]=1;
  		effest[i]=0;
  	}

  	/* Equality of haplotypic effets*/
  	if (hypoth>0)
  	{
  		for (i=0;i<hypoth;i++)
   		{
   			v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
    			v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
    			itp[v2-1]=0;
			if (v1!=1)
			{
				effest[v2-1]=effest[v1-1];
			}
			else
			{
				effest[v2-1]=0;
			}

   		}
  	}
  	/* Homogeneity of allelic effect according to haplotypic background */
  	if (interor==1)
  	{
  		for (j=0;j<nbhypor;j++)
   		{
   			for (i=1;i<nbor[j];i++)
    			{
    				itp[tabinter[j][i][1]-1]=0;
	 			if ( (tabinter[j][0][0]==1) && (tabinter[j][i][0]>1))
	 			{
	 				effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1];
	 			}
	 			else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]>1))
     				{
     					effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]-effest[tabinter[j][0][0]-1];
	 			}
	 			else if ( (tabinter[j][i][0]==1) && (tabinter[j][0][0]==1))
     				{
     					effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1];
	 			}
	 			else
	 			{
	 				effest[tabinter[j][i][1]-1]=effest[tabinter[j][0][1]-1]+effest[tabinter[j][i][0]-1]-effest[tabinter[j][0][0]-1];
	 			}
 			}
   		}
  	}
  	/* Contrainst on haplotype * covariate interactions */
  	if (hypint>0)
  	{
  		for (i=0;i<hypint;i++)
   		{
   			itp[nbhest+ajust+nbadd+tabhypint[i][1]-1]=0;
    			effest[nbhest+ajust+nbadd+tabhypint[i][1]-1]=effest[nbhest+ajust+nbadd+tabhypint[i][0]-1];
   		}
  	}

  	if ((chxt>0) && (haplozero==1))
  	{
  		itp[0]=1;
   		for (i=1;i<nbhest;i++)
   		{
   			itp[i]=0;
   		}
  	}

	j=0;
  	for (i=0;i<maxn;i++)
  	{
  		nitp[i]=-1;
  	}
  	for (i=0;i<nall;i++)
  	{
  		if (itp[i]==1)
  		{
  			j++;
  			nitp[i]=j-1;
  		}
  	}
  	n=j;
  	if (hypoth>0)
  	{
  		for (i=0;i<hypoth;i++)
   		{
   			v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
    			v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
    			if (v1!=1)
    			{
    				nitp[v2-1]=nitp[v1-1];
    			}
			else
			{
				nitp[v2-1]=-2;
			}
   		}
  	}

  	if (hypint>0)
  	{
  		for (i=0;i<hypint;i++)
   		{
   			nitp[nbhest+ajust+nbadd+tabhypint[i][1]-1]=nitp[nbhest+ajust+nbadd+tabhypint[i][0]-1];
   		}
  	}

 	if (chxt==1)
 	{
 		if (msdata==1)
 		{
 			effest[0]=0.5*(log(nbcas)-log(nbtot-nbcas));
 		}
  		else if (msdata==0)
  		{
  			effest[0]=0.5*(log(nbcasm)-log(tabmq[0]-nbcasm));
  		}
 	}
 	else if (chxt==2)
 	{
 		effest[0]=0.5*mean;
 	}
  	if (parfile != NULL) fclose(parfile);parfile=NULL;
        if (readParam != NULL) fclose(readParam);readParam=NULL;

}

/*************************Lecture des contraintes sur les effets pour polytomous*********************/
void lecteffe5()
{
	char rep;
 	int i,j,v1,v2,idc,old;double vv,valf;
 	short chgt;
 	FILE *parfile;
	FILE *readParam;
	div_t biz;
 	haplozero=1;

 	readParam = fopen("paramData.thi","r");

 	printf("Do you want to estimate haplotypic effect (y/n) ?\n");
 	do
 	{
 		fscanf(readParam,"%c",&rep);
 	}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

 	if ((rep=='y') || (rep=='Y'))
 	{
 		haplozero=0;
 	}

 	if ((haplozero==0) && (chxt>0))
 	{
 		hypoth=0;
  		printf("Do you want to test specific hypothesis on haplotypic effects (y/n) ?\n");
  		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

  		if ((rep=='y') || (rep=='Y'))
  		{
  			printf("How many constraints of equality do you wish to consider ?\n");
   			fscanf(readParam,"%d",&hypoth);
   			tabhypo=(int **)malloc((size_t) (hypoth*sizeof(int *)));

   			for (i=0;i<hypoth;i++)
   			{
   				tabhypo[i]=(int *)malloc((size_t) (2*sizeof(int)));
   			}
   			for (i=0;i<hypoth;i++)
   			{
   				printf("Hypothesis %d: enter the two haplotypes ID whose effect must set to be equal\n",i+1);
    				fscanf(readParam,"%d %d",&v1,&v2);
    				tabhypo[i][0]=v1;
    				tabhypo[i][1]=v2;
				if ( (tabhypo[i][0]<0) || (tabhypo[i][1]<0))
				{
					printf("Please check this ID again and Enter some news values\n");
     					i-=1;
    				}
   			}
  		}
  		interor=0;
  		printf("Do you want to test for the homogeneity of some allelic effects  (y/n) ?\n");
  		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));
  		if ((rep=='y') || (rep=='Y'))
  		{
  			interor=1;
   			printf("How many different hypotheses do you wish to test ?\n");
   			do
   			{
   				fscanf(readParam,"%d",&nbhypor);
   			}while (nbhypor>maxhypor);

   			tabinter=(int ***)malloc((size_t) (nbhypor*sizeof(int **)));
   			nbor=(int *)malloc((size_t) (nbhypor*sizeof(int)));
   			for (j=0;j<nbhypor;j++)
   			{
   				chgt=0;
    				printf("Hypothesis %d: How many ORs or Differences to be tested ?\n",j+1);
				do
				{
					fscanf(readParam,"%d",&v1);
					nbor[j]=v1;
				}while (nbor[j]<2);

				tabinter[j]=(int **)malloc((size_t) (nbor[j]*sizeof(int *)));
				for (i=0;i<nbor[j];i++)
				{
					tabinter[j][i]=(int *)malloc((size_t) (2*sizeof(int)));
				}
				for (i=0;i<nbor[j];i++)
				{
					printf("Enter the number of the two haplotypes ID involved in OR/Difference %d\n",i+1);
     					fscanf(readParam,"%d %d",&v1,&v2);
	 				tabinter[j][i][0]=v1;
	 				tabinter[j][i][1]=v2;

	 				if (div(v2,nkat).quot==0)
	 				{
	 					chgt=1;
	 				}

	 				if (chgt==1)
	 				{
	 					for (i=0;i<nbor[j];i++)
	  					{
	  						v1=tabinter[j][i][0];
	   						tabinter[j][i][0]=tabinter[j][i][1];
	   						tabinter[j][i][1]=v1;
      						}
	 				}
	 				if ( (v1<0) || (v2<0))
	 				{
	 					printf("Please check this ID again and Enter some news values\n");
      						i-=1;
     					}
	 				if ((div(v1,nkat).rem) !=(div(v2,nkat).rem))
	 				{
	 					printf("Error in the ID numbers\n");
	 					exit(0);
	 				}
    				}
   			}
  		}
  		nbadd=0;
  		printf("Do you want to test the non-additivity of some haplotypic effects (y/n)\n");
 		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

  		if ((rep=='y') || (rep=='Y'))
  		{
  			printf("How many hypothesis to be tested ?\n");
  			do
  			{
  				fscanf(readParam,"%d",&nbadd);
  			}while (nbadd>maxadd);

  			tadd=(int **)malloc((size_t) (nbadd*sizeof(int *)));
  			for (j=0;j<nbadd;j++)
  			{
  				tadd[j]=(int *)malloc((size_t) (2*sizeof(int)));
   				printf("Hypothesis %d: Enter the two haplotypic effects assumed to interact \n",j+1);
   				fscanf(readParam,"%d %d",&v1,&v2);
   				tadd[j][0]=v1;tadd[j][1]=v2;
  				//  div_t biz; biz=div(fabs(v1-v2),nkat);
  			}
 		}
	}
 	if (ajust>0)
 	{
 		printf("Do you want to test haplotype * environment interactions (y/n) ?");
  		do
  		{
  			fscanf(readParam,"%c",&rep);
  		}while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));

  		if ((rep=='y') || (rep=='Y'))
  		{
  			printf("How many interactions to be tested ?\n");
   			do
   			{
   				fscanf(readParam,"%d",&intercov);
   			}while (intercov>maxcov);

   			tabint=(int **)malloc((size_t) (intercov*sizeof(int *)));

   			for (j=0;j<intercov;j++)
   			{
   				tabint[j]=(int *)malloc((size_t) (2*sizeof(int)));
    				printf("Interaction %d: Enter first the haplotype ID and then the covariate number.\n",j+1);
    				fscanf(readParam,"%d %d",&v1,&v2);
    				tabint[j][0]=v1;
    				tabint[j][1]=v2;
    				if (tabint[j][0]<nkat)
    				{
    					printf("Please check for this Haplotype ID that should be greater than %d\n",nkat-1);
    					exit(0);
    				}
   			}
  		}
 	}

parfile=NULL;
 if ((parfile=fopen("para.txt","r"))==NULL)
 {printf("Error in Para file reading....\n");exit(0);}
  fscanf(parfile,"%d\n",&nbhest);
  nall=(nbhest+ajust+nbadd)*nkat+intercov;


  numhap=(int *)malloc((size_t) (nbhest*sizeof(int)));
  for (i=0;i<nbhhypo;i++) freqest[i]=0;
  for (i=0;i<nbhest;i++)
  {fscanf(parfile,"%d %lf %d %lf\n",&v2,&vv,&idc,&valf);
   for (j=0;j<nkat;j++) {itp[nkat*i+j]=idc;effest[nkat*i+j]=valf;}
   numhap[i]=fcoda2[v2];
   freqest[numhap[i]]=vv;
  }

  vv=0;for (i=0;i<nbhhypo;i++) vv+=freqest[i];
  for (i=0;i<nbhhypo;i++) freqest[i]/=vv;
  vv=0;for (i=0;i<nbhhypo;i++) vv+=freqest[i];

  for (i=nbhest*nkat;i<nall;i++)  {itp[i]=1;effest[i]=0;}

  if (hypoth>0)
  {for (i=0;i<hypoth;i++)
   {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
    v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
    itp[v2]=0;
    if (div(v1,nkat).quot!=0) {effest[v2]=effest[v1];}
                         else {effest[v2]=0;}
   }
  }


  if (haplozero==1)
  {for (j=0;j<nkat;j++) itp[j]=1;
   for (i=nkat;i<nbhest*nkat;i++) itp[i]=0;
  }

  j=0;
  for (i=0;i<nall;i++) {nitp[i]=-1;}
  for (i=0;i<nall;i++) {if (itp[i]==1) {j++;nitp[i]=j-1;}}
  n=j;


  if (hypoth>0)
  {for (i=0;i<hypoth;i++)
   {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
    v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
    if (div(v1,nkat).quot!=0) {nitp[v2]=nitp[v1];}
	else {nitp[v2]=-2;}

   }
  }


  if (interor==1)
  {for (j=0;j<nbhypor;j++)
   {for (i=1;i<nbor[j];i++)
    {itp[tabinter[j][i][1]]=0;
	 if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][i][0],nkat).quot>0))
	 {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]]+effest[tabinter[j][i][0]];
     }
	 else if ( (div(tabinter[j][i][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
     {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]]-effest[tabinter[j][0][0]];
     }
	 else if ( (div(tabinter[j][i][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
     {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]];
     }
	 else
	 {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]]+effest[tabinter[j][i][0]]-effest[tabinter[j][0][0]];
	 }
	}
   }
  }

  if (parfile != NULL) fclose(parfile);

  /*12-10-06 Modif David pour ajouter dans Para.txt le nombre de classe nkat*/
  parfile=NULL;
  if ((parfile=fopen("para.txt","a+"))==NULL)
 {printf("Error in Para file reading....\n");exit(0);}
  printf("nkat: %d\n",nkat);
  fprintf(parfile,"%d\n",nkat);

 if (parfile != NULL) fclose(parfile);parfile=NULL;
 if (readParam != NULL) fclose(readParam);readParam=NULL;
//Fin des modif



}

/*************************Generation aléatoire des haplotypes ambigus*********************/
void generhap()
{double probacum,aleat,quot,dum;
 int i,j,nh=0,nf=0,nh0=0,nh1=0,nf0=0,nf1=0,nls=0;
 short fini;

 for (i=0;i<nbhhypo;i++) {tempfreq[i]=0;}
 if ((chxt==1) || (chxt==3) ||  (chxt==4) ||  (chxt==6)) {for (i=0;i<nbhhypo;i++) {tempdist[i][0]=0;tempdist[i][1]=0;}}

 if (chxt==5) {for (i=0;i<nbhhypo;i++) for (j=0;j<nbcatego;j++) {tempdist[i][j]=0;}}

 if (xlnk==0)
 {suiv=base;
  while ((suiv!=NULL) && (suiv->next!=NULL))
  {if (suiv!=NULL)
   { nls++;
   if (suiv->tnbhapo==1)
   {suiv->hapest[0]=suiv->idnb[0][0];
    suiv->hapest[1]=suiv->idnb[0][1];
   }
   else if (suiv->tnbhapo>1)
   {probacum=0;
    aleat=rand()/(1.0*RAND_MAX);
    fini=0;
    quot=probatot();
    for (i=0;i<suiv->tnbhapo;i++)
    {if (fini==0)
     {dum=probacond(i)/quot;
	  probacum=probacum+dum;
	  if (probacum>aleat)
	  {fini=1;
       suiv->hapest[0]=suiv->idnb[i][0];
       suiv->hapest[1]=suiv->idnb[i][1];
	  }
	 }
    }
   }
   if (suiv->tnbhapo>0)
   {tempfreq[suiv->hapest[0]]+=1;tempfreq[suiv->hapest[1]]+=1;
    if ( (chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
    {tempdist[suiv->hapest[0]][0]+=1*(suiv->phen[0]==0);
     tempdist[suiv->hapest[1]][0]+=1*(suiv->phen[0]==0);
     tempdist[suiv->hapest[0]][1]+=1*(suiv->phen[0]==1);
     tempdist[suiv->hapest[1]][1]+=1*(suiv->phen[0]==1);
    }
	if  (chxt==5)
    {for (i=0;i<nbcatego;i++)
	 {tempdist[suiv->hapest[0]][i]+=1*(suiv->phen[0]==(i+1));
      tempdist[suiv->hapest[1]][i]+=1*(suiv->phen[0]==(i+1));
     }
    }
   }
   suiv=suiv->next;
   }
  }
  for (i=0;i<nbhhypo;i++) {freqest[i]=tempfreq[i]/(2*nbused);}
  if ( (chxt==1) || (chxt==3) || (chxt==4) ||  (chxt==6))
  for (i=0;i<nbhhypo;i++)
  {freqdist[i][0]=tempdist[i][0]/(2*nbtem);freqdist[i][1]=tempdist[i][1]/(2*nbcas);}
  if (chxt==5)
  for (i=0;i<nbhhypo;i++) {for (j=0;j<nbcatego;j++) freqdist[i][j]=tempdist[i][j]/(2*nbsujktgo[j]);}
 }
 else if (xlnk==1)
 {suiv=base;
  nh=0;nf=0;
  while ((suiv!=NULL) && (suiv->next!=NULL))
  {if (((int) suiv->z[0])==0)
   {
	if (suiv->tnbhapo==1)
    {nh+=1;
     suiv->hapest[0]=suiv->idnb[0][0];
     suiv->hapest[1]=suiv->idnb[0][1];
     tempfreq[suiv->hapest[0]]+=1;
     if ( (chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
     {tempdist[suiv->hapest[0]][0]+=1*(suiv->phen[0]==0);
      nh0+=(suiv->phen[0]==0);
      tempdist[suiv->hapest[0]][1]+=1*(suiv->phen[0]==1);
	  nh1+=(suiv->phen[0]==1);
     }
	}
    //else {printf("Option not handled \n");exit(0);}
   }
   else if (((int) suiv->z[0])==1)
   {if (suiv->tnbhapo==1)
    {suiv->hapest[0]=suiv->idnb[0][0];
     suiv->hapest[1]=suiv->idnb[0][1];
    }
    else if (suiv->tnbhapo>1)
    {probacum=0;
     aleat=rand()/(1.0*RAND_MAX);
     fini=0;
     quot=Xprobatot();
     for (i=0;i<suiv->tnbhapo;i++)
     {if (fini==0)
      {dum=Xprobacond(i)/quot;
	   probacum=probacum+dum;
	   if (probacum>aleat)
	   {fini=1;
        suiv->hapest[0]=suiv->idnb[i][0];
        suiv->hapest[1]=suiv->idnb[i][1];
	   }
	  }
     }
    }
    if (suiv->tnbhapo>0)
    {tempfreq[suiv->hapest[0]]+=1;tempfreq[suiv->hapest[1]]+=1;
     nf+=1;
     if ( (chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
     {tempdist[suiv->hapest[0]][0]+=1*(suiv->phen[0]==0);
      tempdist[suiv->hapest[1]][0]+=1*(suiv->phen[0]==0);
      tempdist[suiv->hapest[0]][1]+=1*(suiv->phen[0]==1);
      tempdist[suiv->hapest[1]][1]+=1*(suiv->phen[0]==1);
      nf0+=(suiv->phen[0]==0);
	  nf1+=(suiv->phen[0]==1);
     }
	 if (chxt==5) {printf("PARTIE A REMPLIR \n");exit(0);}
    }
   }
   suiv=suiv->next;
  }
  for (i=0;i<nbhhypo;i++) {freqest[i]=tempfreq[i]/(2*nf+nh);}
  if ( (chxt==1) || (chxt==3) || (chxt==4) ||  (chxt==6))
  for (i=0;i<nbhhypo;i++)
  {freqdist[i][0]=tempdist[i][0]/(2*nf0+nh0);freqdist[i][1]=tempdist[i][1]/(2*nf1+nh1);}
  if (chxt==5) {printf("PARTIE A REMPLIR \n");exit(0);}
 }


 suiv=NULL;
}

double llambda(double temps)
{double res,res1,val1;
 int h1,hh1,h2,hh2,i,j;
 dindividu *pv1,*pv2;

 pv1 = base;
 res1=0;res=0;
 while ((pv1!=NULL) && (pv1->next!=NULL))
 {if ( (pv1->tnbhapo>0) && (pv1->phen[0]==1.0) && (pv1->phen[1]<=temps))
  {pv2=base;
   res1=0;
   while ((pv2!= NULL) && (pv2->next!=NULL))
   {if ((pv2->phen[1]>=pv1->phen[1]) && (pv2->tnbhapo>0))
	{val1=0;
	 for (i=0;i<ajust;i++) val1+=effest[nbhest+i]*pv2->z[i];
     if (haplozero==0)
     {h1=pv2->hapest[0];h2=pv2->hapest[1];
	  hh1=coding(h1);   hh2=coding(h2);
      if (hh1>0) {val1+=effest[hh1];}
      if (hh2>0) {val1+=effest[hh2];}
      if (nbadd>0)
      {for (j=0;j<nbadd;j++)
       {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
        {val1+=effest[nbhest+ajust+j];}
       }
      }
      for (j=0;j<intercov;j++)
      {val1+=pv2->z[tabint[j][1]-1]*effest[nbhest+ajust+nbadd+j]*((hh1==tabint[j][0]-1)+(hh2==tabint[j][0]-1));
      }
	 }
	 res1+=exp(val1);
    }
	pv2=pv2->next;
   } /*pv2*/
   res+=(1/res1);
  }/*if*/
  pv1=pv1->next;
 }/*pv1*/

// New DAVID
   free ((dindividu *) pv1);free ((dindividu *) pv2);pv1=NULL;pv2=NULL;
 return (res);
}


/**************************Calcul des probabilités conditionnelles*****************************/
double probacond(int val)
{double res=1.0,val1,val2,vv2,likeli,ss=0,ls=0,ps=0,num,denom;
 int h1,h2,hh1,hh2,idx,i,j;


 idx=1*(haplozero==1)+nbhest*(haplozero==0);
 h1=suiv->idnb[val][0];
 h2=suiv->idnb[val][1];

 if (chxt==0) {res*=freqest[h1]*freqest[h2]*(2-(h1==h2));}
 else if (chxt<5)
 { val1=suiv->phen[0];val2=2*effest[0];
   for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
   if ((chxt==1) && (offset==1)) val2+=suiv->z[ajust];
   if (haplozero==0)
   {hh1=coding(h1);hh2=coding(h2);
     if (hh1>0) {val2+=effest[hh1];}
     if (hh2>0) {val2+=effest[hh2];}
     if (nbadd>0)
     {for (i=0;i<nbadd;i++)
      {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
       {val2+=effest[nbhest+ajust+i];}
      }
     }
     for (i=0;i<intercov;i++)
     {val2+=suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
     }
    }
    if ((chxt==1) || (chxt==4))
    {res*=freqest[h1]*freqest[h2]*(2-(h1==h2))*exp(val2*val1)/(1+exp(val2));
    }
    else if (chxt==2)
   {vv2=-0.5*(val1-val2)*(val1-val2)/(ste*ste);
    likeli=exp(vv2)/(ste*sqrt(2.0*pi));
    res*=freqest[h1]*freqest[h2]*(2-(h1==h2))*likeli;
   }
   else if (chxt==3)
   {res*=freqest[h1]*freqest[h2]*(2-(h1==h2));
     if (suiv->phen[0]==0)
     {res*=exp(-exp(val2));
     }
     else if (suiv->phen[0]==1)
     {res*=exp(val2)*exp(-exp(val2));
     }
    }
  /*else if (chxt==6)
  {res*=freqest[h1]*freqest[h2]*(2-(h1==h2));
   ss=exp(-exp(val2)*pow(effest[nall-1]*suiv->phen[1],effest[nall-2]));
   ls=effest[nall-1]*effest[nall-2]*pow(suiv->phen[1]*effest[nall-1],effest[nall-2]-1)*exp(val2);
   ps=ss*ls;
   res*=ss;
   if (suiv->phen[0]==1) {res*=ls;}
    correction for P(S=1)
   res/=(suiv->wgt+(1-suiv->wgt)*ps);
  }*/
 }
 else if (chxt==5)
 {idx= ((int) suiv->phen[0])-1;
  val2=0;

  denom=1;
  for (i=0;i<nkat;i++)
  {val2=2*effest[i];
   for (j=0;j<ajust;j++) val2+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];
   if (haplozero==0)
   {hh1=coding(h1);hh2=coding(h2);
    if (hh1>0) {val2+=effest[hh1*nkat+i];}
    if (hh2>0) {val2+=effest[hh2*nkat+i];}
	for (j=0;j<nbadd;j++)
    {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
      {val2+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
    }
	for (j=0;j<intercov;j++)
    {val2+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
    }

   }
   denom+=exp(val2);
  }
  if (idx==0) num=1;
  else if (idx>0)
  {val2=2*effest[idx-1];
   for (j=0;j<ajust;j++) val2+=effest[nbhest*nkat+j*nkat+idx-1]*suiv->z[j];
   if (haplozero==0)
   {hh1=coding(h1);hh2=coding(h2);
    if (hh1>0) {val2+=effest[nkat*hh1+idx-1];}
    if (hh2>0) {val2+=effest[nkat*hh2+idx-1];}
	for (j=0;j<nbadd;j++)
    {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
     {val2+=effest[(nbhest+ajust)*nkat+j*nkat+idx-1];}
    }
	for (j=0;j<intercov;j++)
    {val2+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+idx-1)==tabint[j][0])+ ((hh2*nkat+idx-1)==tabint[j][0]));
    }
   }
   num=exp(val2);
  }
  res*=freqest[h1]*freqest[h2]*(2-(h1==h2))*num/denom;
 }
 return(res);
}



/**************************Calcul des probabilités conditionnelles*****************************/
double Xprobacond(int val)
{double res=1.0,val1,val2,vv2,likeli,ss=0,ls=0,ps=0;
 int h1,h2,hh1,hh2,idx,i;
 idx=1*(haplozero==1)+nbhest*(haplozero==0);

 h1=suiv->idnb[val][0];
 h2=suiv->idnb[val][1];

 if (chxt==0) {res*=freqest[h1]*freqest[h2]*(2-(h1==h2));}
 else
 {val1=suiv->phen[0];val2=effest[0];
  for (i=0;i<ajust;i++) val2+=effest[nbhest+i]*suiv->z[i];
  if ((chxt==1) && (offset==1)) val2+=suiv->z[ajust];
  if (haplozero==0)
  {hh1=coding(h1);hh2=coding(h2);
   if (hh1>0) {val2+=0.5*effest[hh1];}
   if (hh2>0) {val2+=0.5*effest[hh2];}
   if (nbadd>0)
   {for (i=0;i<nbadd;i++)
    {if ( ((tadd[i][0]-1==hh1) && (tadd[i][1]-1==hh2)) || ((tadd[i][0]-1==hh2) && (tadd[i][1]-1==hh1)) )
     {val2+=effest[nbhest+ajust+i];}
    }
   }
   for (i=0;i<intercov;i++)
   {val2+=0.5*suiv->z[tabint[i][1]-1]*effest[nbhest+ajust+nbadd+i]*((hh1==tabint[i][0]-1)+(hh2==tabint[i][0]-1));
   }
  }
  if ((chxt==1) || (chxt==4))
  {res*=freqest[h1]*freqest[h2]*(2-(h1==h2))*exp(val2*val1)/(1+exp(val2));
  }
  else if (chxt==2)
  {vv2=-0.5*(val1-val2)*(val1-val2)/(ste*ste);
   likeli=exp(vv2)/(ste*sqrt(2.0*pi));
   res*=freqest[h1]*freqest[h2]*(2-(h1==h2))*likeli;
  }
  else if (chxt==3)
  {res*=freqest[h1]*freqest[h2]*(2-(h1==h2));
   if (suiv->phen[0]==0)
   {res*=exp(-exp(val2));
   }
   else if (suiv->phen[0]==1)
   {res*=exp(val2)*exp(-exp(val2));
   }
  }
  else if (chxt==6)
  {res*=freqest[h1]*freqest[h2]*(2-(h1==h2));
   ss=exp(-exp(val2)*pow(effest[nall-1]*suiv->phen[1],effest[nall-2]));
   ls=effest[nall-1]*effest[nall-2]*pow(suiv->phen[1]*effest[nall-1],effest[nall-2]-1)*exp(val2);
   ps=ss*ls;
   res*=ss;
   if (suiv->phen[0]==1) {res*=ls;}
   /* correction for P(S=1)*/
   res/=(suiv->wgt+(1-suiv->wgt)*ps);
  }
 }
 return(res);
}




/**************************Calcul des probabilités totales************************************/
double probatot()
{double som;
 int j;
 som=0;
 for (j=0;j<(suiv->tnbhapo);j++)
 {som+=probacond(j);}
 return(som);
}

double Xprobatot()
{double som;
 int j;
 som=0;
 for (j=0;j<suiv->tnbhapo;j++)
 {som+=Xprobacond(j);}
 return(som);
}






/**************************Programme Principal*********************************************/
//int main()
int thesiasRun(char *fileName,
	       int pMaxvarfic,
	       int pNbloci,
	       int *pIdloci,
	       int pLdmatrix,
	       int pMsdata,
	       int pR2,
	       int pChxt,
	       int pNum0,
	       int pIdtime,
	       int pOffset,
    	       int pIdoffset,
    	       int pAjust,
    	       int *pNumajust,
    	       int pXlnk,
	       int pNumsx)
{int i,mynit,j,k,ii,ll,imn,imnk,idn,idk,ll2;
 char l,rep;
 short trouve,ik,il;
 int maxhap,*numord,*numobs,*place,v1,v2,nbwald,nc, noUse;
 double *freqcas,*freqtem,somt,res,*frqobs,*frqord,vraise,vraiscond,valse,valef;
 double **ldp,**lr2,**chidp,mink,minl,dmax;
 double *hwc,*bresl;
 double **vinter,**vintra;
 double **freqkat;
 int *idwald;
 div_t biz;

 double timedif;
 double time1 = (double)clock();

 matrixp matse;
 time_t now;
 FILE *outres,*outfile;
//VG 14112006 : initialisation des variables globales
 interor=0;hypoth=0;nbhypor=0;nbadd=0;hypint=0;intercov=0;xlnk=0;
 maxhapair=0;nbloci=0; nbhhypo=0;nnt=0;ajust=0;nbcatego=0;nkat=0;nall=0;n=0;
 
 
 //lecture();
 printf("DEBUT\n");
 noUse = lecture(fileName,
	     	 pMaxvarfic,
	     	 pNbloci,
	     	 pIdloci,
	     	 pLdmatrix,
	     	 pMsdata,
	     	 pR2,
	     	 pChxt,
	     	 pNum0,
	     	 pIdtime,
	     	 pOffset,
	     	 pIdoffset,
	     	 pAjust,
	     	 pNumajust,
	     	 pXlnk,
	     	 pNumsx);


 maxhap=ipow(2,nbloci);
 fcoda1=(short *) malloc((size_t) (maxhap*sizeof(short)));
 fcoda2=(int *) malloc((size_t) (maxhap*sizeof(int)));
 for (i=0;i<maxhap;i++) {fcoda1[i]=0;fcoda2[i]=-1;}


 effest=(double *) malloc ((size_t) (maxn*sizeof(double)));
 itp=(int *) malloc ((size_t) (maxn*sizeof(int)));
 nitp=(int *) malloc ((size_t) (maxn*sizeof(int)));
 itptp=(int *) malloc ((size_t) (maxn*sizeof(int)));
 nitptp=(int *) malloc ((size_t) (maxn*sizeof(int)));
 hwc=(double *) malloc ((size_t) (nbloci*sizeof(double)));

printf("FIN LECTURE\n");
 determhapo();
 printf("FIN FCT1\n");
 hapopres();

 distrmq();


 if (chxt==5) {categorie();}


 initfreq(hwc);
 recodage();

 if ((chxt>0) && (chxt!=5)) {lecteffe();}
 if (chxt==5) {lecteffe5();}


 modif2=(double *) malloc((size_t) ((nall)*sizeof(double )));
 for (i=0;i<nall;i++) modif2[i]=0;
 mdvs2=(double *) malloc((size_t) ((nall)*sizeof(double )));
 for (i=0;i<nall;i++) mdvs2[i]=0;
 dmat2=(double *) malloc((size_t) ((nall)*sizeof(double )));
 for (i=0;i<nall;i++) dmat2[i]=0;


 if (chxt==3)
 {tritime();
  tabpi=(double *) malloc((size_t) ((nbtotused)*sizeof(double )));
  for (i=0;i<nbtotused;i++) tabpi[i]=0;
  tablo=(double **) malloc((size_t) ((nbtotused)*sizeof(double *)));
  for (i=0;i<nbtotused;i++) tablo[i]=(double *) malloc((size_t) ((nall)*sizeof(double)));
  for (i=0;i<nbtotused;i++) {for (j=0;j<nall;j++) tablo[i][j]=0;}

  //;initlist();plus necessaire

 }
 if (chxt==4)
 {tripair();}



 moyfreq=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
 if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
 {freqcas=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
  freqtem=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
 }
 if (chxt==5)
 {freqkat=(double **) malloc((size_t) (nbcatego*sizeof(double *)));
  for (i=0;i<nbcatego;i++)
  {freqkat[i]=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
   for (j=0;j<nbhhypo;j++) freqkat[i][j]=0.0;
  }
 }


 inclus=(short *) malloc((size_t) (nbhhypo*sizeof(short)));

 for (i=0;i<nbhhypo;i++) {moyfreq[i]=0;inclus[i]=0;}
 if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
 {for (i=0;i<nbhhypo;i++) {freqcas[i]=0;freqtem[i]=0;}
 }

 printf("RUNNING THE SEM ALGORITHM ....\n");
 srand(time(&now)%40);


 moyeff=(double *) malloc ((size_t) (nall*sizeof(double)));
 for (i=0;i<nall;i++) {moyeff[i]=0;}

 if ( (chxt==6) || (chxt==3) || (chxt==5) || (chxt==4))
 {vecbeta=(double **) malloc ((size_t) (nall*sizeof(double *)));
  for (i=0;i<nall;i++) vecbeta[i]=(double *) malloc ((size_t) ((mymaxit-nburn-1)*sizeof(double)));
  for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbeta[i][j]=0;
  vecbetat=(double **) malloc ((size_t) ((mymaxit-nburn-1)*sizeof(double *)));
  for (i=0;i<mymaxit-nburn-1;i++) vecbetat[i]=(double *) malloc ((size_t) (nall*sizeof(double)));
  for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbetat[j][i]=0;

  vintra=(double **) malloc ((size_t) (nall*sizeof(double *)));
  for (i=0;i<nall;i++) vintra[i]=(double *) malloc ((size_t) (nall*sizeof(double )));
  for (i=0;i<nall;i++) for (j=0;j<nall;j++) vintra[i][j]=0;

  tabres=(double *)malloc((size_t) (3*sizeof(double)));
  for (i=0;i<3;i++) tabres[i]=0;
 }



 if (chxt==0)
 {for (mynit=1;mynit<mymaxit;mynit++)
  {generhap();
   presence();
   if (mynit%100==0) {printf("Iteration %d\n", mynit);}
   if (mynit==nburn) {printf("The burning period has finished....\n");}
   if (mynit>nburn)
   {for (i=0;i<nbhhypo;i++) {moyfreq[i]+=freqest[i];}
   }
  }
 }
 else if (chxt!=5)
 {for (mynit=1;mynit<mymaxit;mynit++)
  {generhap();
   presence();
   for (i=0;i<nall;i++) {itptp[i]=itp[i];nitptp[i]=nitp[i];}
   j=0;
   for (i=0;i<nbhest;i++)
   {k=numhap[i];
    if (inclus[k]==0) {itptp[i]=0;nitptp[i]=-1;}
    if ( (chxt==1) && (inclus[k]==1) && ( (freqdist[k][0]==0) || (freqdist[k][1]==0)))
    {itptp[i]=0;nitptp[i]=-1;}
   }
   j=0;
   for (i=0;i<nall;i++) {if (itptp[i]==1) {j++;nitptp[i]=j-1;}}
   n=j;
   if (hypoth>0)
   {for (i=0;i<hypoth;i++)
    {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
     v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
     if (v1!=1) {nitptp[v2-1]=nitptp[v1-1];}
     else {nitptp[v2-1]=-2;}
    }
   }
   if (hypint>0)
   {for (i=0;i<hypint;i++)
    {nitptp[nbhest+ajust+nbadd+tabhypint[i][1]-1]=nitptp[nbhest+ajust+nbadd+tabhypint[i][0]-1];}
   }
   for (i=0;i<nall;i++) {effest[i]=0;}
   conti=go;
   while (conti !=quit)
   {tempx=0;
    if (chxt<3)  {if (xlnk==1) Xfisherscoring();
	              else  fisherscoring();
				 }
    if (chxt==3) {/*coxtempo();*/coxtablo();}
    if (chxt==4) {matchpair();}

    if (tempx<seuilfisher) {conti=quit;}
   }
   if (mynit%100==0) {printf("Iteration %d\n", mynit);}
   if (mynit==nburn) {printf("The burning period has finished....\n");}
   if (mynit>nburn)
   {if (chxt==2) meanste+=ste;
 	else if (chxt==3)
    {for (i=0;i<nall;i++) vecbeta[i][mynit-nburn-1]=effest[i];
     coxrubin(vintra);
	}
    else if (chxt==4)
	{/* Formulation LOUIS*/
     for (i=0;i<nall;i++) vecbeta[i][mynit-nburn-1]=effest[i];
  	 pairedfish(vintra);

	}
    for (i=0;i<nbhhypo;i++) {moyfreq[i]+=freqest[i];}
    if  (chxt!=2) for (i=0;i<nbhhypo;i++) {freqtem[i]+=freqdist[i][0];freqcas[i]+=freqdist[i][1];}
    for (i=0;i<nall;i++) moyeff[i]+=effest[i];
   }

  }
 }
 else if (chxt==5)
 {for (mynit=1;mynit<mymaxit;mynit++)
  {generhap();
   presence();
   for (i=0;i<nall;i++) {itptp[i]=itp[i];nitptp[i]=nitp[i];}
   j=0;

   for (i=0;i<nbhest;i++)
   {k=numhap[i];
    if (inclus[k]==0) for (j=0;j<nkat;j++) {itptp[i*nkat+j]=0;nitptp[i*nkat+j]=-1;}
	if (inclus[k]==1)
    {ii=0;
     /*ATTENTION A CETTE VERIFICATION , IL FAUT QUE L'HAPLOTYPE  EXISTE DANS LA POPULA de REFE?*/
	 for (j=0;j<nkat;j++)
	 {if ( (freqdist[k][j+1]==0) || (freqdist[k][0]==0))
      {itptp[i*nkat+j]=0;nitptp[i*nkat+j]=-1;
	  }
     }
    }
   }
   j=0;
   for (i=0;i<nall;i++) {if (itptp[i]==1) {j++;nitptp[i]=j-1;}}
   n=j;

   if (hypoth>0)
   {for (i=0;i<hypoth;i++)
    {v1=tabhypo[i][0]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][1]*(tabhypo[i][1]<tabhypo[i][0]);
     v2=tabhypo[i][1]*(tabhypo[i][0]<tabhypo[i][1])+tabhypo[i][0]*(tabhypo[i][1]<tabhypo[i][0]);
     if (div(v1,nkat).quot!=0) {nitptp[v2]=nitptp[v1];}
     else {nitptp[v2]=-2;}
    }
   }
   /*
   if (hypint>0)
   {for (i=0;i<hypint;i++)
    {nitptp[nbhest+ajust+nbadd+tabhypint[i][1]-1]=nitptp[nbhest+ajust+nbadd+tabhypint[i][0]-1];}
   }*/

   for (i=0;i<nall;i++) {effest[i]=0;}
   conti=go;
   while (conti !=quit)
   {tempx=0;
    polytomous();
	//if (xlnk==1) Xfisherscoring();  else  fisherscoring();
   	if (tempx<seuilfisher) {conti=quit;}
   }
   if (mynit%100==0) {printf("Iteration %d\n", mynit);}
   if (mynit==nburn) {printf("The burning period has finished....\n");}
   if (mynit>nburn)
   {for (i=0;i<nall;i++) vecbeta[i][mynit-nburn-1]=effest[i];
    vpolyto(vintra);
   	for (i=0;i<nbhhypo;i++) {moyfreq[i]+=freqest[i];}
    for (i=0;i<nbhhypo;i++) for (j=0;j<nbcatego;j++) {freqkat[j][i]+=freqdist[i][j];}
	for (i=0;i<nall;i++) moyeff[i]+=effest[i];
   }
  }
 }

 if (chxt==2) {meanste/=(mymaxit-nburn-1);}
 for (i=0;i<nall;i++) {moyeff[i]/=(mymaxit-nburn-1);}
 for (i=0;i<nbhhypo;i++) {moyfreq[i]/=(mymaxit-nburn-1);}
 if ((chxt==1) || (chxt==3) || (chxt==4) || (chxt==6))
	 for (i=0;i<nbhhypo;i++) {freqtem[i]/=(mymaxit-nburn-1);freqcas[i]/=(mymaxit-nburn-1);}
 if (chxt==5) for (i=0;i<nbhhypo;i++) for (j=0;j<nbcatego;j++) {freqkat[j][i]/=(mymaxit-nburn-1);}




 n=0;
 for (i=0;i<nall;i++) {effest[i]=moyeff[i];  if (itp[i]==1) n++;  }

 outfile=NULL;
 if ((outfile=fopen("result.htm","w"))==NULL)
 {printf("Error in Output file writing....\n");exit(0);}

 outres=NULL;
 if ((outres=fopen("result.txt","w"))==NULL)
 {printf("Error in Output file writing....\n");exit(0);}




 fprintf(outres,"Total number of individuals: %d\n\n",nbtot);

// fprintf(outfile,"<BODY BGCOLOR=""#CCFFCC"">\n");
 fprintf(outfile,"<BODY BGCOLOR=""#F9E6C4"">\n");
 fprintf(outfile,"<table align=center border=0 width=80%%>\n");
 fprintf(outfile,"<tr><td width = 50%% align=left>Total number of individuals:</td>");
 fprintf(outfile,"<td align=left>%d</td></tr>\n",nbtot);

 for (i=0;i<nbloci+1;i++)
 {fprintf(outfile,"<tr><td align=left width=50%%>Number of individuals with %d missing data:</td>\n",i);
  fprintf(outfile,"<td align=left >%d</td></tr>\n",tabmq[i]);
  fprintf(outres,"Number of individuals with %d missing data: %d\n",i,tabmq[i]);
 }

 if (xlnk==1)
 {fprintf(outres,"\nTotal number of males/females: %d / %d\n",nbhf[0][0],nbhf[1][0]);
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=left width=50%%>Total number of males/females </td>\n");
  fprintf(outfile,"<td>%d / %d </td></tr>\n",nbhf[0][0],nbhf[1][0]);
 }

 fprintf(outres,"\n");
 fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");

 if ((chxt==1) || (chxt==4))
 {fprintf(outfile,"<tr><td align=left width=50%%>Total number of cas </td>\n");
  fprintf(outfile,"<td>%d</td></tr>\n",nbcas);
  fprintf(outfile,"<tr><td align=left width=50%%>Total number of cas without missing data </td>\n");
  fprintf(outfile,"<td>%d</td></tr>\n",nbcasm);
  fprintf(outres,"Total number of cas: %d\n",nbcas);
  fprintf(outres,"Total number of cas without missing data: %d\n",nbcasm);


  if (xlnk==1)
  {fprintf(outres,"\nTotal number of males /females in cases: %d / %d\n",nbhf[0][2],nbhf[1][2]);
   fprintf(outres,"\nTotal number of males /females in controls: %d / %d\n",nbhf[0][1],nbhf[1][1]);

   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=left width=50%%>Total number of males/females in cases</td>\n");
   fprintf(outfile,"<td>%d / %d </td></tr>\n",nbhf[0][2],nbhf[1][2]);
   fprintf(outfile,"<tr><td align=left width=50%%>Total number of males/females in controls</td>\n");
   fprintf(outfile,"<td>%d / %d </td></tr>\n",nbhf[0][1],nbhf[1][1]);
  }





 }
 else if (chxt==2)
 {fprintf(outfile,"<tr><td align=left width=50%%>Phenotypic Mean = %f</td>\n",mean);
  fprintf(outfile,"<td align=left>Standard Error  = %f</td></tr>\n\n",ste0);
  fprintf(outres,"Phenotypic Mean = %f\n",mean);
  fprintf(outres,"Standard Error  = %f\n\n",ste0);
 }
 else if ( (chxt==3) || (chxt==6))
 {fprintf(outfile,"<tr><td align=left width=50%%>Total number of uncensored individuals </td>\n");
  fprintf(outfile,"<td>%d</td></tr>\n",nbcas);
  fprintf(outfile,"<tr><td align=left width=50%%>Total number of uncensored individuals without missing data </td>\n");
  fprintf(outfile,"<td>%d</td></tr>\n",nbcasm);
  fprintf(outres,"Total number of uncensored individuals: %d\n",nbcas);
  fprintf(outres,"Total number of uncensored individuals without missing data: %d\n",nbcasm);
 }
 else if (chxt==5)
 {fprintf(outfile,"<tr>");
  for (i=0;i<nbcatego;i++)
  {fprintf(outfile,"<td align=left width=30%%>Phenotypic Level %d</td>",i+1);
   fprintf(outfile,"<td align=left width=70%%>%d</td>",nbsujktgo[i]);
   fprintf(outfile,"</tr>");
  }
  fprintf(outres,"\n");
  for (i=0;i<nbcatego;i++)
  {fprintf(outres,"Phenotypic Level %d\t%d",i+1,nbsujktgo[i]);
   fprintf(outres,"\n");
  }
 }


 fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr></table>\n");
 fprintf(outres,"\n\n");


 fprintf(outfile,"<br><table align=center border=0 width=80%%>\n");
 for (i=0;i<nbloci;i++)
 {fprintf(outfile,"<tr><td align=left width=40%%>Allele frequency at locus %d </td>\n",i+1);
  fprintf(outfile,"<td align=left>(%c/%c)</td>\n",letter[i][0],letter[i][1]);
  fprintf(outres,"Allele frequency at locus %d ",i+1);
  fprintf(outres,"(%c/%c) ",letter[i][0],letter[i][1]);
  fprintf(outfile,"<td align=left colspan=""2"">%.5f / %.5f   p(HWE) = %f</td></tr>\n",alfreq[i],1-alfreq[i],chdtrc(1.0,hwc[i]));
  fprintf(outres,"%.5f / %.5f\tp(HWE) = %f\n",alfreq[i],1-alfreq[i],chdtrc(1.0,hwc[i]));

 }

 fprintf(outres,"\n");
 fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr></table>\n");

 fprintf(outres,"Frequencies of Plausible Haplotype under Linkage Equilibrium\n\n");

 fprintf(outfile,"<table align=center border=0 width=80%%>\n");
 fprintf(outfile,"<tr>\n<td align=center colspan=""5"">Frequencies of Plausible Haplotype under Linkage Equilibrium</td></tr>\n\n");
 fprintf(outfile,"<tr><td> </td></tr>\n");

 vect1=tnbhbase;i=0;
 while (vect1!=NULL)
 {if (vect1->present==1)
  {fprintf(outres,"Haplotype [%d]\t",i);
   fprintf(outfile,"<tr><td align=left width=20%%>Haplotype [%d]</td>\n",i);
   fprintf(outfile,"<td align=center width =25%%>");
   for (j=0;j<nbloci;j++)
   {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
    fprintf(outfile,"%c",l);
    fprintf(outres,"%c",l);
   }
   fprintf(outfile,"</td><td align=left colspan=""3"">%f</td></tr>\n",vect1->frqle);
   fprintf(outres,"\t%f\n",vect1->frqle);

  }
  i++;
  vect1=vect1->down;
 }
 fprintf(outres,"\n\n");
 fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");



  if (chxt==2) {ste=meanste;}
  if (chxt==4)
  {bresl=(double *)malloc((size_t) (3*sizeof(double)));
   likematchpair(moyfreq,moyeff,bresl);
   nnt=nbhest+n-1;
/*  fishpair(moyfreq,moyeff,matse); ANCIENNE FORMULATION*/
   vinter=(double **) malloc ((size_t) (nall*sizeof(double *)));
   for (i=0;i<nall;i++) vinter[i]=(double *) malloc ((size_t) ((nall)*sizeof(double)));
   for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbeta[i][j]-=effest[i];
   for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbetat[j][i]=vecbeta[i][j];
   for (i=0;i<nall;i++)
	for (j=0;j<nall;j++)
	   {vinter[i][j]=0;
        for (k=0;k<mymaxit-nburn-1;k++)
		  vinter[i][j]+=(vecbeta[i][k])*(vecbetat[k][j]);
	   }
   for (i=0;i<nall;i++) for (j=0;j<nall;j++) vinter[i][j]*=(mymaxit-nburn)/((mymaxit-nburn-1)*(mymaxit-nburn-2));


   for (i=0;i<nall;i++) for (j=0;j<nall;j++) vinter[i][j]+=vintra[i][j]/(mymaxit-nburn-1);
   for (i=0;i<nall;i++)
 	 for (j=0;j<nall;j++) matse[nbhest-1+nitp[i]-1][nbhest-1+nitp[j]-1]=vinter[i][j];
 }
  else if ( (chxt==6) || (chxt==3))
  {if (chxt==3)
   {bresl=(double *)malloc((size_t) (3*sizeof(double)));
    breslow1(moyfreq,moyeff,bresl);
    nnt=nbhest+n-1;
   }
   vinter=(double **) malloc ((size_t) (nall*sizeof(double *)));
   for (i=0;i<nall;i++) vinter[i]=(double *) malloc ((size_t) ((nall)*sizeof(double)));
   for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbeta[i][j]-=effest[i];
   for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbetat[j][i]=vecbeta[i][j];
   for (i=0;i<nall;i++)
	for (j=0;j<nall;j++)
	   {vinter[i][j]=0;
        for (k=0;k<mymaxit-nburn-1;k++)
		  vinter[i][j]+=(vecbeta[i][k])*(vecbetat[k][j]);
	   }
  for (i=0;i<nall;i++) for (j=0;j<nall;j++) vinter[i][j]*=(mymaxit-nburn)/((mymaxit-nburn-1)*(mymaxit-nburn-2));
  for (i=0;i<nall;i++) for (j=0;j<nall;j++) vinter[i][j]+=vintra[i][j]/(mymaxit-nburn-1);
  for (i=0;i<nall;i++)
	 for (j=0;j<nall;j++) matse[nbhest-1+nitp[i]-1][nbhest-1+nitp[j]-1]=vinter[i][j];
  }
  else if (chxt==5)
  {vinter=(double **) malloc ((size_t) (nall*sizeof(double *)));
   for (i=0;i<nall;i++) vinter[i]=(double *) malloc ((size_t) ((nall)*sizeof(double)));
   for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbeta[i][j]-=effest[i];
   for (i=0;i<nall;i++) for (j=0;j<mymaxit-nburn-1;j++) vecbetat[j][i]=vecbeta[i][j];
   for (i=0;i<nall;i++)
	for (j=0;j<nall;j++)
	   {vinter[i][j]=0;
        for (k=0;k<mymaxit-nburn-1;k++)
		  vinter[i][j]+=(vecbeta[i][k])*(vecbetat[k][j]);
	   }
   for (i=0;i<nall;i++) for (j=0;j<nall;j++) vinter[i][j]*=(mymaxit-nburn)/((mymaxit-nburn-1)*(mymaxit-nburn-2));
   for (i=0;i<nall;i++) for (j=0;j<nall;j++) vinter[i][j]+=vintra[i][j]/(mymaxit-nburn-1);
   for (i=0;i<nall;i++)
	 for (j=0;j<nall;j++) matse[nbhest-1+nitp[i]][nbhest-1+nitp[j]]=vinter[i][j];

   vraise=likepoly(moyfreq,moyeff);
   vraiscond=vraise-condlike(moyfreq);
  }
  else if (chxt>0) {if (xlnk==1) {Xfishem(moyfreq,moyeff,matse);}
	                else fishem(moyfreq,moyeff,matse);
				   }
  else if (chxt==0) {if (xlnk==1) {Xfishnull(moyfreq,matse);}
	                 else fishnull(moyfreq,matse);
					}

  place=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
  for (i=0;i<nbhhypo;i++) place[i]=-1;
  j=0;for (i=0;i<nbhhypo;i++) if (moyfreq[i]>0) {place[i]=j;j++;}

 vect1=NULL;
 fprintf(outfile,"<tr>\n<td align=center colspan=""5"">After %d iterations burn at step %d</td></tr>\n\n",mymaxit,nburn);
 fprintf(outfile,"<tr><td> </td></tr>\n");
 fprintf(outres,"After %d iterations burn at step %d\n\n",mymaxit,nburn);


 if (chxt==0)
 {fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype Frequencies under Linkage Disequilibrium</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=left> </td><td align=left></td><td align=left>Estimation</td><td align=left>StError</td>\n");
  fprintf(outfile,"<td align=left>T-Test</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");

  fprintf(outres,"Frequencies of Plausible Haplotypes under Linkage Disequilibrium\n\n");
  fprintf(outres,"\t\t\t\tEstimation\tStError\tT-Test\n\n");

  vect1=tnbhbase;i=0;somt=0;
  while (vect1!=NULL)
  {ll=fcoda2[vect1->numnew];
   if (vect1->present==1)
   {if (moyfreq[ll]>0)
    {fprintf(outfile,"<tr><td align=left width=20%%>Haplotype [%d] </td>\n",i);
     fprintf(outfile,"<td align=center width =25%%>");
     fprintf(outres,"Haplotype [%d] \t",i);
	 for (j=0;j<nbloci;j++)
     {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
      fprintf(outfile,"%c",l);
      fprintf(outres,"%c",l);
     }
	 fprintf(outfile,"</td><td align=left>");
     res=moyfreq[ll];
     fprintf(outfile,"%f</td>",res);
	 fprintf(outres," %f\t",res);
	 if ((place[ll]!=0) && (res>0))
     {valse=sqrt(matse[place[ll]-1][place[ll]-1]);
      fprintf(outfile,"<td align=left>%f</td><td align=left>%f</td></tr>\n",valse,res/valse);
      fprintf(outres,"%f\t%f\n",valse,res/valse);
     }
     else {fprintf(outfile,"</tr>\n");
	       fprintf(outres,"\n");
	   }
    }
   }
   i++;
   vect1=vect1->down;
  }
  vect1=NULL;
 }

 if (chxt==0)
 {frqobs=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
  frqord=(double *) malloc((size_t) (nbhhypo*sizeof(double)));
  numobs=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
  numord=(int *) malloc((size_t) (nbhhypo*sizeof(int)));
  for (i=0;i<nbhhypo;i++) {frqobs[i]=0;frqord[i]=0;numord[i]=0;numobs[i]=0;}
  j=0;
  vect1=tnbhbase;
  while (vect1!=NULL)
  {if (vect1->present==1)
   {i=fcoda2[vect1->numnew];
    if (moyfreq[i]>0) {frqobs[j]=moyfreq[i]/*/(mymaxit-nburn-1)*/;numobs[j]=vect1->numnew;j++;}
   }
   vect1=vect1->down;
  }
  for (i=0;i<j;i++)
  {if (i==0) {numord[i]=numobs[i];frqord[i]=frqobs[i];}
   else
   {trouve=0;k=0;
    while ((trouve==0) && (k<i))
	{if (frqobs[i]>frqord[k]) {trouve=1;}
     k++;
    }
  	if (trouve==0) {numord[i]=numobs[i];frqord[i]=frqobs[i];}
    else if (trouve==1)
	{for (ll=i;ll>k-1;ll--) {numord[ll]=numord[ll-1];frqord[ll]=frqord[ll-1];}
	 numord[k-1]=numobs[i];frqord[k-1]=frqobs[i];
    }
   }
  }
   parafile=NULL;
  if ((parafile=fopen("para.txt","w"))==NULL)
  {printf("Error in Output file writing....\n");exit(0);}
  fprintf(parafile,"%d\n",j);
  for (i=0;i<j;i++)
  {if (i==0) fprintf(parafile,"%d %f 1 %f\n",numord[i],frqord[i],effest[0]);
   else fprintf(parafile,"%d %f 0 0\n",numord[i],frqord[i]);
  }


  if (parafile != NULL) fclose(parafile);

 }



 if ((chxt==1) || (chxt==2) )
 {fprintf(outres,"\n");
  fprintf(outres,"Estimated Haplotype Frequencies under Linkage Disequilibrium\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotypes Frequencies under Linkage Disequilibrium</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td> </td><td align=left>Estimation</td><td align=left>StError</td>\n");
  fprintf(outfile,"<td align=left>T-Test</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outres,"\n");
  fprintf(outres,"\t\t\t\tEstimation\tStError\t\tT-Test\n\n");
  vect1=tnbhbase;i=0;somt=0;
  while (vect1!=NULL)
  {ll=fcoda2[vect1->numnew];
   if (vect1->present==1)
   {if (moyfreq[ll]>0)
    {fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d [%d]</td>\n",coding(ll)+1,i);
     fprintf(outres,"Haplotype %d [%d] \t",coding(ll)+1,i);
     fprintf(outfile,"<td align=center width =25%%>");
     for (j=0;j<nbloci;j++)
     {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
      fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
     }
     fprintf(outfile,"</td><td align=left>");
     res=moyfreq[ll];
     fprintf(outfile,"\t%f",res);fprintf(outres,"\t%f",res);
     if ((place[ll]!=0) && (res>0))
     {valse=sqrt(matse[place[ll]-1][place[ll]-1]);
      fprintf(outfile,"</td><td align=left>%f</td><td align=left>%f</td></tr>\n",valse,res/valse);
      fprintf(outres,"\t%f\t%f\n",valse,res/valse);
     }
     else {fprintf(outfile,"</td></tr>\n");
           fprintf(outres,"\n");
          }
     somt+=res;
    }
   }
   i++;
   vect1=vect1->down;
  }
  vect1=NULL;
 }
 else if ((chxt==3) || (chxt==4) || (chxt==6) || (chxt==5))
 {fprintf(outres,"\n");
  fprintf(outres,"Estimated Haplotype Frequencies under Linkage Disequilibrium\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotypes Frequencies under Linkage Disequilibrium</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td> </td><td align=left>Estimation</td><td align=left>&nbsp</td>\n");
  fprintf(outfile,"<td align=left>&nbsp </td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outres,"\n");
  fprintf(outres,"\t\t\t\tEstimation\t\t\t\n\n");
  vect1=tnbhbase;i=0;somt=0;
  while (vect1!=NULL)
  {ll=fcoda2[vect1->numnew];
   if (vect1->present==1)
   {if (moyfreq[ll]>0)
    {fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d [%d]</td>\n",coding(ll)+1,i);
     fprintf(outres,"Haplotype %d [%d] \t",coding(ll)+1,i);
     fprintf(outfile,"<td align=center width =25%%>");
     for (j=0;j<nbloci;j++)
     {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
      fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
     }
     fprintf(outfile,"</td><td align=left>");
     res=moyfreq[ll];
     fprintf(outfile,"\t%f",res);fprintf(outres,"\t%f",res);
     fprintf(outfile,"</td></tr>\n");
     fprintf(outres,"\n");
     somt+=res;
    }
   }
   i++;
   vect1=vect1->down;
  }
  vect1=NULL;
 }

 if ((chxt==1))
 {fprintf(outres,"\n\n");
  fprintf(outres,"Estimated Haplotype Frequencies in\tControls (=0)\tCases (=1)\n\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td  align=left>Haplotype Frequencies in</td><td align=left>Controls (=0)</td>\n");
  fprintf(outfile,"<td  align=left>Cases (=1)</td><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  vect1=tnbhbase;i=0;
  while (vect1!=NULL)
  {ll=fcoda2[vect1->numnew];
   if (vect1->present==1)
   {if (moyfreq[ll]>0)
	{fprintf(outres,"Haplotype %d [%d] \t",coding(ll)+1,i);
	 fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d [%d]</td>\n",coding(ll)+1,i);
     fprintf(outfile,"<td align=center width =25%%>");
	 for (j=0;j<nbloci;j++)
     {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
      fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
     }
	 fprintf(outfile,"</td><td align=left>");

     fprintf(outfile,"%f</td><td align=left>%f</td><td> </td></tr>\n",freqtem[ll],freqcas[ll]);
     fprintf(outres,"\t\t%f\t%f\n",freqtem[ll],freqcas[ll]);
	}
   }
   i++;
   vect1=vect1->down;
  }
  vect1=NULL;
 }
 else if ((chxt==3)  || (chxt==4) || (chxt==6) )
 {fprintf(outres,"\n\n");
  fprintf(outres,"Estimated Haplotype Frequencies in\tCensored (=0)\tUncensored (=1) individuals\n\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td  align=left>Haplotype Frequencies in</td><td align=left>Censored (=0)</td>\n");
  fprintf(outfile,"<td  align=left>Uncensored (=1)</td><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  vect1=tnbhbase;i=0;
  while (vect1!=NULL)
  {ll=fcoda2[vect1->numnew];
   if (vect1->present==1)
   {if (moyfreq[ll]>0)
	{fprintf(outres,"Haplotype %d [%d] \t",coding(ll)+1,i);
	 fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d [%d]</td>\n",coding(ll)+1,i);
     fprintf(outfile,"<td align=center width =25%%>");
	 for (j=0;j<nbloci;j++)
     {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
      fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
     }
	 fprintf(outfile,"</td><td align=left>");
     fprintf(outfile,"%f</td><td align=left>%f</td><td> </td></tr>\n",freqtem[ll],freqcas[ll]);
     fprintf(outres,"\t\t%f\t%f\n",freqtem[ll],freqcas[ll]);
	}
   }
   i++;
   vect1=vect1->down;
  }
  vect1=NULL;
 }
 else if (chxt==5)
 {fprintf(outres,"\n\n");
  fprintf(outres,"Estimated Haplotype Frequencies in\n\n\t\t\t\t");
  for (i=0;i<nkat+1;i++) fprintf(outres,"Class level %d\t",i);fprintf(outres,"\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td  colspan=""4"" align=center> Haplotype Frequencies </td></tr>");
  fprintf(outfile,"<tr><td>&nbsp</td><td>&nbsp</td>\n");
  for (i=0;i<nkat+1;i++) fprintf(outfile,"<td align=left>Class Level %d</td>",i+1);
  fprintf(outfile,"</tr>\n");
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");

  vect1=tnbhbase;i=0;
  while (vect1!=NULL)
  {ll=fcoda2[vect1->numnew];
   if (vect1->present==1)
   {if (moyfreq[ll]>0)
	{fprintf(outres,"Haplotype %d [%d] \t",coding(ll)+1,i);
	 fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d [%d]</td>\n",coding(ll)+1,i);
     fprintf(outfile,"<td align=center width =25%%>");
	 for (j=0;j<nbloci;j++)
     {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
      fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
     }
	 fprintf(outfile,"</td>");
     for (j=0;j<nbcatego;j++)
	 {fprintf(outfile,"<td align=left>%lf</td>",freqkat[j][ll]);
      fprintf(outres,"\t%f",freqkat[j][ll]);
	 }
	 fprintf(outfile,"</tr>");fprintf(outres,"\n");

	}
   }
   i++;
   vect1=vect1->down;
  }
  vect1=NULL;
 }


 fprintf(outres,"\n\n");
 fprintf(outfile,"<tr><td> </td></tr>\n");

 if ((chxt==1) || (chxt==2))
 {fprintf(outres,"\n");
  fprintf(outres,"\t\t\tHaplotype Effects\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype Effects</td></tr>\n");
  if ((chxt==1) && (haplozero==0))
  {fprintf(outfile,"<tr><td align=center colspan=""5"">(Haplotypic OR by comparison to the reference with its 95%% CI)</td></tr>\n\n");
   fprintf(outres,"\t\t(Haplotypic OR by comparison to the reference with its 95%% CI)\n\n");
  }
  if ((chxt==2) && (haplozero==0))
  {fprintf(outfile,"<tr><td align=center colspan=""5"">(Haplotypic effect by comparison to the reference with its 95%% CI)</td></tr>\n");
   fprintf(outres,"\t\t(Haplotypic effect by comparison to the reference with its 95%% CI)\n\n");
  }
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td> </td><td align=left>Estimation</td><td align=left>StError</td>\n");
  fprintf(outfile,"<td align=left>T-Test</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outres,"\t\t\t\tEstimation\tStError\t\tT-Test\n\n");

  if (haplozero==1)
  {valse=sqrt(matse[nbhest-1][nbhest-1]);
   fprintf(outfile,"<tr><td align=left>Intercept</td><td> </td><td align=left>%f</td>\n",effest[0]);
   fprintf(outfile,"<td align=left>%f</td><td> </td></tr>\n",valse);
   fprintf(outres,"Intercept\t%f\t%f\n",effest[0],valse);
  }
  else if (haplozero==0)
  {vect1=tnbhbase;trouve=0;
   j=0;k=0;
   while (vect1!=NULL)
   {trouve=0;
    i=0;
	while ((trouve==0) && (i<nbhest))
    {if (fcoda2[vect1->numnew]==numhap[i])
     {
	  trouve=1;
	  if (i==0)
	  {fprintf(outfile,"<tr><td align=left width=20%%>Intercept  1 </td>\n");
       fprintf(outres,"Intercept 1\t");
      }
	  else
	  {fprintf(outres,"Haplotype %d\t",i+1);
	   fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d </td>",i+1);
	  }
	  fprintf(outfile,"<td align=center width =25%%>");
	  for (j=0;j<nbloci;j++)
      {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
       fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
	  }

      fprintf(outfile,"</td>");fprintf(outres,"\t");
      if (i==0) {valse=sqrt(matse[nbhest-1][nbhest-1]);

				 fprintf(outfile,"<td align=left>%.5f</td><td align=left>%f</td><td align=left>%f</td></tr>\n\n",effest[0],valse,effest[0]/valse);
				 fprintf(outfile,"<tr><td> </td></tr>\n");
				 fprintf(outres,"\t%f\t%f\n\n",effest[0],valse,effest[0]/valse);
				}
	  else      {
	            fprintf(outfile,"<td align=left>%.5f</td>",effest[i]);
 	             fprintf(outres,"\t%.5f\t",effest[i]);
			     if (itp[i]==1)
				 {valse=sqrt(matse[nbhest-1+nitp[i]][nbhest-1+nitp[i]]);
			      valef=effest[i];
			      affichage(outfile,valef,valse);
				  affichage2(outres,valef,valse);
				 }
				 else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
				       fprintf(outres,"\n");
				      }
			    }


	 }
     i++;
	}
    vect1=vect1->down;trouve=0;
   }
  }
 }
 else if ( (chxt==3) || (chxt==4) || (chxt==6) )
 {fprintf(outres,"\n");
  fprintf(outres,"\t\t\tHaplotype Effects\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype Effects</td></tr>\n");
  if ((chxt==1) && (haplozero==0))
  {fprintf(outfile,"<tr><td align=center colspan=""5"">(Haplotypic OR by comparison to the reference with its 95%% CI)</td></tr>\n\n");
   fprintf(outres,"\t\t(Haplotypic OR by comparison to the reference with its 95%% CI)\n\n");
  }
  if ((chxt==2) && (haplozero==0))
  {fprintf(outfile,"<tr><td align=center colspan=""5"">(Haplotypic effect by comparison to the reference with its 95%% CI)</td></tr>\n");
   fprintf(outres,"\t\t(Haplotypic effect by comparison to the reference with its 95%% CI)\n\n");
  }
  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td> </td><td align=left>Estimation</td><td align=left>StError</td>\n");
  fprintf(outfile,"<td align=left>T-Test</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outres,"\t\t\t\tEstimation\tStError\t\tT-Test\n\n");

  if (haplozero==1)
  {
  }
  else if (haplozero==0)
  {vect1=tnbhbase;trouve=0;
   j=0;k=0;
   while (vect1!=NULL)
   {trouve=0;
    i=0;
	while ((trouve==0) && (i<nbhest))
    {if (fcoda2[vect1->numnew]==numhap[i])
     {
	  trouve=1;
	  if (i==0)
	  {fprintf(outfile,"<tr><td align=left width=20%%>Intercept  1 </td>\n");
       fprintf(outres,"Intercept 1\t");
      }
	  else
	  {fprintf(outres,"Haplotype %d\t",i+1);
	   fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d </td>",i+1);
	  }
	  fprintf(outfile,"<td align=center width =25%%>");
	  for (j=0;j<nbloci;j++)
      {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
       fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
	  }

      fprintf(outfile,"</td>");fprintf(outres,"\t");
      if (i==0) {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
				 fprintf(outres,"\n");
				}
	  else      {
	            fprintf(outfile,"<td align=left>%.5f</td>",effest[i]);
 	             fprintf(outres,"\t%.5f\t",effest[i]);
                 if (itp[i]==1)
				 {valse=sqrt(matse[nbhest-1+nitp[i]-1][nbhest-1+nitp[i]-1]);
			      valef=effest[i];
			      affichage(outfile,valef,valse);
				  affichage2(outres,valef,valse);
				 }
				 else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
				       fprintf(outres,"\n");
				      }
			    }


	 }
     i++;
	}
    vect1=vect1->down;trouve=0;
   }
  }
 }
 else if (chxt==5)
 {fprintf(outres,"\n");
  fprintf(outres,"\t\t\tHaplotype Effects\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype Effects</td></tr>\n");

  fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td><td> </td><td align=left>Estimation</td><td align=left>StError</td>\n");
  fprintf(outfile,"<td align=left>T-Test</td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outres,"\t\t\t\tEstimation\tStError\t\tT-Test\n\n");

  for (nc=0;nc<nkat;nc++)
  {fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outres,"\n\n");
   if (haplozero==1)
   {valse=sqrt(matse[nbhest-1+nc][nbhest-1+nc]);
    fprintf(outfile,"<tr><td align=left width=30%%>Intercept for Class Level %d</td><td> </td><td align=left>%f</td>\n",nc+2,effest[nc]);
    fprintf(outfile,"<td align=left>%f</td><td> </td></tr>\n",valse);
    fprintf(outres,"Intercept for Class Level %d\t%f\t%f\n",nc+2,effest[nc],valse);
   }
   else if (haplozero==0)
   {fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotypic OR for Class Level %d by comparison to Class Level 1</td></tr>\n\n",nc+2);
    fprintf(outres,"\t\tHaplotypic OR for Class Level %d by comparison to Class Level 1\n\n",nc+2);

    vect1=tnbhbase;trouve=0;
    j=0;k=0;
    while (vect1!=NULL)
    {trouve=0;
     i=0;
	 while ((trouve==0) && (i<nbhest))
     {if (fcoda2[vect1->numnew]==numhap[i])
      {trouve=1;
	   if (i==0)
	   {fprintf(outfile,"<tr><td align=left width=20%%>Intercept (ID:%d) </td>\n",nc);
        fprintf(outres,"Intercept (ID:%d)\t",nc);
       }
	   else
	   {fprintf(outres,"Haplotype %d (ID:%d)\t",i+1,i*nkat+nc);
	    fprintf(outfile,"<tr><td align=left width=20%%>Haplotype %d (ID:%d)</td>",i+1,i*nkat+nc);
	   }
	   fprintf(outfile,"<td align=center width =25%%>");
	   for (j=0;j<nbloci;j++)
       {l=letter[j][0]*(vect1->listall[j]==1)+letter[j][1]*(vect1->listall[j]==2);
        fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
	   }
       fprintf(outfile,"</td>");fprintf(outres,"\t");
       if (i==0) {valse=sqrt(matse[nbhest-1+nc][nbhest-1+nc]);
 				  fprintf(outfile,"<td align=left>%.5f</td><td align=left>%f</td><td align=left>%f</td></tr>\n\n",effest[nc],valse,effest[nc]/valse);
				  fprintf(outfile,"<tr><td> </td></tr>\n");
				  fprintf(outres,"\t%f\t%f\n\n",effest[nc],valse,effest[nc]/valse);
				 }
	   else      {fprintf(outfile,"<td align=left>%.5f</td>",effest[i*nkat+nc]);
 	              fprintf(outres,"\t%.5f\t",effest[i*nkat+nc]);
			      if (itp[i*nkat+nc]==1)
				  {valse=sqrt(matse[nbhest-1+nitp[i*nkat+nc]][nbhest-1+nitp[i*nkat+nc]]);
			       valef=effest[i*nkat+nc];
			       affichage(outfile,valef,valse);
				   affichage2(outres,valef,valse);
				  }
				  else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
				        fprintf(outres,"\n");
				       }
			     }


	  }
      i++;
	 }
     vect1=vect1->down;trouve=0;
    }
   }
  }
  fprintf(outres,"\n\n");
  fprintf(outres,"Log-likelihood = %f ",vraise);
  fprintf(outres,"(Nb studied subjects = %d)\n",nbused);
  fprintf(outfile,"<tr><td> </td></tr>\n");fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align =left colspan=""2"">Log-likelihood %f</td>\n",vraise);
  fprintf(outfile,"<td align=left colspan=""3"">(Nb studied subjects = %d)</td></tr>",nbused);
  fprintf(outfile,"<tr><td> </td></tr>\n");

  fprintf(outres,"Conditionnal Log-likelihood = %f ",vraiscond);
  fprintf(outres,"(df = %2d)\n\n\n",n);
  fprintf(outfile,"<tr><td align=left colspan=""2"">Conditionnal Log-likelihood=%f</td>\n",vraiscond);
  fprintf(outfile,"<td align=left colspan=""3"">(df = %2d)</td></tr>",n);
  fprintf(outfile,"<tr><td> </td></tr>\n");

  if (ajust>0)
  {fprintf(outres,"\n\n");
   fprintf(outres,"\t\t\tCovariable Adjustment\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Covariable Adjustment</td></tr>\n");
   for (nc=0;nc<nkat;nc++)
   {fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outfile,"<tr><td colspan=""5"">Information for CLASS LEVEL %d </td></tr>\n",nc+2);
    fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outres,"\n");
    fprintf(outres,"\n");
    fprintf(outres,"Information for CLASS LEVEL %d\n",nc+2);
    fprintf(outres,"\n");
    for (i=0;i<ajust;i++)
    {fprintf(outfile,"<tr><td align=left width=20%%>Covariate %d</td>\n",i+1);
     fprintf(outfile,"<td align=left width=25%%>Column Number %d</td>\n",numajust[i]);
     fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest*nkat+i*nkat+nc]);
     fprintf(outres,"Covariate %d \n",i+1);
     fprintf(outres,"Column Number %d \t\t",numajust[i]);
     fprintf(outres,"%f\t",effest[nbhest*nkat+i*nkat+nc]);
     if (itp[nbhest*nkat+i*nkat+nc]==1)
     {valse=sqrt(matse[nbhest-1+nitp[nbhest*nkat+i*nkat+nc]][nbhest-1+nitp[nbhest*nkat+i*nkat+nc]]);
      valef=effest[nbhest*nkat+i*nkat+nc];
	  affichage(outfile,valef,valse);affichage2(outres,valef,valse);
     }
     else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
           fprintf(outres,"\n");
          }
    }
   }
  }
  if (nbadd>0)
  {fprintf(outres,"\n\n");
   fprintf(outres,"\t\t\tNon additivity of haplotypic effects\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Non additivity of haplotypic effects</td></tr>\n");

   for (nc=0;nc<nkat;nc++)
   {fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outfile,"<tr><td colspan=""5"">Information for CLASS LEVEL %d </td></tr>\n",nc+2);
    fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outres,"\n");
    fprintf(outres,"\n");
    fprintf(outres,"Information for CLASS LEVEL %d\n",nc+2);
    fprintf(outres,"\n");

    for (i=0;i<nbadd;i++)
    {fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>Haplotype %d x Haplotype %d </td>\n",tadd[i][0],tadd[i][1]);
     fprintf(outfile,"<td align=left>%f</td>\n",effest[(nbhest+ajust)*nkat+i*nkat+nc]);
     fprintf(outres,"Haplotype %d x Haplotype %d\t",tadd[i][0],tadd[i][1]);
     fprintf(outres,"%f\t",effest[(nbhest+ajust)*nkat+i*nkat+nc]);
     if (itp[(nbhest+ajust)*nkat+i*nkat+nc]==1)
     {valse=sqrt(matse[nbhest-1+nitp[(nbhest+ajust)*nkat+i*nkat+nc]][nbhest-1+nitp[(nbhest+ajust)*nkat+i*nkat+nc]]);
      valef=effest[(nbhest+ajust)*nkat+i*nkat+nc];
	  affichage(outfile,valef,valse);affichage2(outres,valef,valse);
     }
     else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
           fprintf(outres,"\n");
          }
    }
   }
  }
  if (intercov>0)
  {fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype x Covariable Interaction</td></tr>");
   fprintf(outres,"\n\n");
   fprintf(outres,"\t\t\tHaplotype x Covariable Interaction\n\n");
   for (i=0;i<intercov;i++)
   {biz=div(tabint[i][0],nkat);
    fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>CLASS LEVEL %d Haplotype %d with Covariate %d </td>\n",biz.rem+2,biz.quot+1,tabint[i][1]);
    fprintf(outres,"CLASS LEVEL %d Haplotype %d with Covariate %d\t",biz.rem+2,biz.quot+1,tabint[i][1]);
    fprintf(outfile,"<td align=left>%f</td>\n",effest[(nbhest+ajust+nbadd)*nkat+i]);
    fprintf(outres,"%f\t",effest[(nbhest+ajust+nbadd)*nkat+i]);
    if (itp[(nbhest+ajust+nbadd)*nkat+i]==1)
    {valse=sqrt(matse[nbhest-1+nitp[(nbhest+ajust+nbadd)*nkat+i]][nbhest-1+nitp[(nbhest+ajust+nbadd)*nkat+i]]);
     valef=effest[(nbhest+ajust+nbadd)*nkat+i];
     affichage(outfile,valef,valse);affichage2(outres,valef,valse);
	}
    else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");fprintf(outres,"\n");}
   }
  }
 }



 if (chxt==2)
 {fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
  fprintf(outfile,"<tr><td align=left colspan=""5"">Global Standard Error   %f</td></tr>\n",ste0);
  fprintf(outfile,"<tr><td align=left colspan=""5"">Residual Standard Error %f</td></tr>\n",ste);
  fprintf(outres,"\n\n");
  fprintf(outres,"Global Standard Error   %f\n",ste0);
  fprintf(outres,"Residual Standard Error %f\n",ste);
 }
 if ( (chxt==1) || (chxt==2))
 {if (ajust>0)
  {fprintf(outres,"\n");
   fprintf(outres,"\t\t\tCovariable Adjustment\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Covariable Adjustment</td></tr>\n");
   for (i=0;i<ajust;i++)
   {if ((i==0) && (xlnk==1))
    {fprintf(outfile,"<tr><td align=left width=20%%>Gender </td>\n");
     fprintf(outres,"Gender %d \n",i+1);
    }
    else
    {fprintf(outfile,"<tr><td align=left width=20%%>Covariate %d</td>\n",i+1);
     fprintf(outres,"Covariate %d \n",i+1);
    }
    fprintf(outfile,"<td align=left width=25%%>Column Number %d</td>\n",numajust[i]);
    fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest+i]);
    fprintf(outres,"Column Number %d \t\t",numajust[i]);
    fprintf(outres,"%f\t",effest[nbhest+i]);
    if (itp[nbhest+i]==1)
    {valse=sqrt(matse[nbhest-1+nitp[nbhest+i]][nbhest-1+nitp[nbhest+i]]);
     valef=effest[nbhest+i];
 	 affichage(outfile,valef,valse);affichage2(outres,valef,valse);
    }
    else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
          fprintf(outres,"\n");
         }
   }
   if (intercov>0)
   {fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype x Covariable Interaction</td></tr>");
    fprintf(outres,"\n\n");
    fprintf(outres,"\t\t\tHaplotype x Covariable Interaction\n\n");
    for (i=0;i<intercov;i++)
    {if ((xlnk==1) && (tabint[i][1]==1))
	 {fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>Haplotype %d x Gender </td>\n",tabint[i][0]);
      fprintf(outres,"Haplotype %d x Gender\t",tabint[i][0]);
     }
	 else
	 {fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>Haplotype %d x Covariate %d </td>\n",tabint[i][0],tabint[i][1]);
      fprintf(outres,"Haplotype %d x Covariate %d\t",tabint[i][0],tabint[i][1]);
     }
	 fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest+ajust+i+nbadd]);
     fprintf(outres,"%f\t",effest[nbhest+ajust+i+nbadd]);
     if (itp[nbhest+ajust+nbadd+i]==1)
     {valse=sqrt(matse[nbhest-1+nitp[nbhest+ajust+nbadd+i]][nbhest-1+nitp[nbhest+ajust+nbadd+i]]);
      valef=effest[nbhest+ajust+nbadd+i];
      affichage(outfile,valef,valse);affichage2(outres,valef,valse);
	 }
     else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");fprintf(outres,"\n");}
    }
   }
  }
  if (nbadd>0)
  {fprintf(outres,"\n\n");
   fprintf(outres,"\t\t\tNon additivity of haplotypic effects\n\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Non additivity of haplotypic effects</td></tr>\n");
   for (i=0;i<nbadd;i++)
   {fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>Haplotype %d x Haplotype %d </td>\n",tadd[i][0],tadd[i][1]);
    fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest+ajust+i]);
    fprintf(outres,"Haplotype %d x Haplotype %d\t",tadd[i][0],tadd[i][1]);
    fprintf(outres,"%f\t",effest[nbhest+ajust+i]);
    if (itp[nbhest+ajust+i]==1)
    {valse=sqrt(matse[nbhest-1+nitp[nbhest+ajust+i]][nbhest-1+nitp[nbhest+ajust+i]]);
     valef=effest[nbhest+ajust+i];
     affichage(outfile,valef,valse);affichage2(outres,valef,valse);
    }
    else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");fprintf(outres,"\n");}
   }
  }
 }
 else if ((chxt==3) || (chxt==4) || (chxt==6))
 {if (ajust>0)
  {fprintf(outres,"\n");
   fprintf(outres,"\t\t\tCovariable Adjustment\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Covariable Adjustment</td></tr>\n");
   for (i=0;i<ajust;i++)
   {fprintf(outfile,"<tr><td align=left width=20%%>Covariate %d</td>\n",i+1);
    fprintf(outfile,"<td align=left width=25%%>Column Number %d</td>\n",numajust[i]);
    fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest+i]);
    fprintf(outres,"Covariate %d \n",i+1);
    fprintf(outres,"Column Number %d \t\t",numajust[i]);
    fprintf(outres,"%f\t",effest[nbhest+i]);
    if (itp[nbhest+i]==1)
    {valse=sqrt(matse[nbhest-1+nitp[nbhest+i]-1][nbhest-1+nitp[nbhest+i]-1]);
     valef=effest[nbhest+i];
	 affichage(outfile,valef,valse);affichage2(outres,valef,valse);
    }
    else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");
          fprintf(outres,"\n");
         }
   }
   if (intercov>0)
   {fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    fprintf(outfile,"<tr><td align=center colspan=""5"">Haplotype x Covariable Interaction</td></tr>");
    fprintf(outres,"\n\n");
    fprintf(outres,"\t\t\tHaplotype x Covariable Interaction\n\n");
    for (i=0;i<intercov;i++)
    {fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>Haplotype %d x Covariate %d </td>\n",tabint[i][0],tabint[i][1]);
     fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest+ajust+i+nbadd]);
     fprintf(outres,"Haplotype %d x Covariate %d\t",tabint[i][0],tabint[i][1]);
     fprintf(outres,"%f\t",effest[nbhest+ajust+i+nbadd]);
     if (itp[nbhest+ajust+nbadd+i]==1)
     {valse=sqrt(matse[nbhest-1+nitp[nbhest+ajust+nbadd+i]-1][nbhest-1+nitp[nbhest+ajust+nbadd+i]-1]);
      valef=effest[nbhest+ajust+nbadd+i];
      affichage(outfile,valef,valse);affichage2(outres,valef,valse);
 	 }
     else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");fprintf(outres,"\n");}
    }
   }
  }
  if (nbadd>0)
  {fprintf(outres,"\n\n");
   fprintf(outres,"\t\t\tNon additivity of haplotypic effects\n\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""5"">Non additivity of haplotypic effects</td></tr>\n");
   for (i=0;i<nbadd;i++)
   {fprintf(outfile,"<tr><td align=left colspan=""2"" width=45%%>Haplotype %d x Haplotype %d </td>\n",tadd[i][0],tadd[i][1]);
    fprintf(outfile,"<td align=left>%f</td>\n",effest[nbhest+ajust+i]);
    fprintf(outres,"Haplotype %d x Haplotype %d\t",tadd[i][0],tadd[i][1]);
    fprintf(outres,"%f\t",effest[nbhest+ajust+i]);
    if (itp[nbhest+ajust+i]==1)
    {valse=sqrt(matse[nbhest-1+nitp[nbhest+ajust+i]-1][nbhest-1+nitp[nbhest+ajust+i]-1]);
     valef=effest[nbhest+ajust+i];
     affichage(outfile,valef,valse);affichage2(outres,valef,valse);
    }
    else {fprintf(outfile,"<td> </td><td> </td><td> </td></tr>\n");fprintf(outres,"\n");}
   }
  }
 }
 if (chxt>0)
 {if (chxt<3)
  {if (xlnk==1)
   {vraise=Xlikelihood(moyfreq,moyeff);
    vraiscond=vraise-Xcondlike(moyfreq);
   }
   else if (xlnk==0)
   {vraise=likelihood(moyfreq,moyeff);
    vraiscond=vraise-condlike(moyfreq);
   }
   fprintf(outres,"\n\n");
   fprintf(outres,"Log-likelihood = %f ",vraise);
   fprintf(outres,"(Nb studied subjects = %d)\n",nbused);
   fprintf(outres,"Conditionnal Log-likelihood = %f ",vraiscond);
   fprintf(outres,"(df = %2d)\n\n\n",nnt-nbhest+1);
   fprintf(outfile,"<tr><td> </td></tr>\n");fprintf(outfile,"<tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align =left colspan=""2"">Log-likelihood %f</td>\n",vraise);
   fprintf(outfile,"<td align=left colspan=""3"">(Nb studied subjects = %d)</td></tr>",nbused);
   fprintf(outfile,"<tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=left colspan=""2"">Conditionnal Log-likelihood=%f</td>\n",vraiscond);
   fprintf(outfile,"<td align=left colspan=""3"">(df = %2d)</td></tr>",nnt-nbhest+1);
   fprintf(outfile,"<tr><td> </td></tr></table>\n");
  }
  else if ( (chxt==3) || (chxt==4))
  {fprintf(outres,"\n\n");
   fprintf(outres,"-2 x Log-likelihood (with covariates)    = %f ",2*bresl[2]);
   fprintf(outres,"(Nb studied subjects = %d)\n",nbused);
   fprintf(outres,"-2 x Log-likelihood (without covariates) = %f ",2*bresl[0]);
   fprintf(outres,"    (df = %2d)\n\n\n",nnt-nbhest+1-1);
   fprintf(outfile,"<tr><td> </td></tr>\n");fprintf(outfile,"<tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align =left colspan=""3"">-2x Log-likelihood (with covariates) = %f</td>\n",2*bresl[2]);
   fprintf(outfile,"<td align=left colspan=""2"">(Nb studied subjects = %d)</td></tr>",nbused);
   fprintf(outfile,"<tr><td> </td></tr>\n");
   fprintf(outfile,"<tr><td align=left colspan=""3"">-2 x Log-likelihood (without covariates) =%f</td>\n",2*bresl[0]);
   fprintf(outfile,"<td align=left colspan=""2"">(df = %2d)</td></tr>",nnt-nbhest+1-1);
   fprintf(outfile,"<tr><td> </td></tr></table>\n");
  }
  if ((haplozero==0) && (chxt!=5))
  {fprintf(outfile,"<table border=0  width=80%%>\n");
   for (i=1;i<nbloci+1;i++)
   {fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outfile,"<tr><td align=left colspan=""5"">Polymorphism %d %c/%c</td></tr>\n",i,letter[i-1][0],letter[i-1][1]);
    fprintf(outres,"\n");
    fprintf(outres,"Polymorphism %d %c/%c\n",i,letter[i-1][0],letter[i-1][1]);
	k=pow(2,nbloci-i);
    j=0;
    vect1=tnbhbase;
    while ( (vect1!=NULL) )
    {ll=fcoda2[vect1->numnew];
     ll2=fcoda2[vect1->numnew+k];
	 if ((vect1->present==1) && (ll2>-1) && (vect1->listall[i-1]==1) && (coding(ll)>-1) && (coding(ll2)>-1) && (moyfreq[ll]>0.00) &&  (moyfreq[ll2]!=0.00))
     {fprintf(outfile,"<tr><td colspan=""2"" align=center>Haplotypic Background</td>\n");
      fprintf(outfile,"<td align=left>");
	  fprintf(outres,"Haplotypic Background ");
	  for (ii=0;ii<nbloci;ii++)
      {if (ii+1!=i)
       {l=letter[ii][0]*(vect1->listall[ii]==1)+letter[ii][1]*(vect1->listall[ii]==2);
        fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
	   }
       else if (ii+1==i) {fprintf(outfile,"-");
	                      fprintf(outres,"-");
					      }
      }

	  fprintf(outfile,"</td><td align= left width=10%%>");
	  fprintf(outfile," %d-%d \t",coding(fcoda2[vect1->numnew])+1,coding(fcoda2[vect1->numnew+k])+1);
	  fprintf(outfile,"</td><td colspan=""2"">");
	  fprintf(outres,"\t%d-%d\t",coding(fcoda2[vect1->numnew])+1,coding(fcoda2[vect1->numnew+k])+1);

	  imn=coding(fcoda2[vect1->numnew]); imnk=coding(fcoda2[vect1->numnew+k]);

	  if ((imn==0) && (itp[imnk]==1))
	  {valef=effest[imnk];
       if ( (chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]]);}
	   else  {valse=sqrt(matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1]);}
      }
      else if ((imn==0) && (itp[imnk]==0))
	  {valef=effest[imnk];
       if (valef!=0)
	   {if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]]);}
	    else   {valse=sqrt(matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1]);}
       }
       else valse=-1;
      }
	  else if ((imnk==0) && (itp[imn]==1))
	  {valef=-effest[imn];
       if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]]);}
	   else   {valse=sqrt(matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1]);}
      }
      else  if ((itp[imn]==1) && (itp[imnk]==1))
	  {if ((imn!=0) && (imnk!=0))
       {if ((chxt!=3) && (chxt!=6) && (chxt!=4))
		{valef=matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]];
	     valef+=matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]];
	     valef-=2*matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imnk]];
        }
		else
		{valef=matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1];
	     valef+=matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1];
	     valef-=2*matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imnk]-1];
        }
		valse=sqrt(valef);
        valef=effest[imnk]-effest[imn];
	   }
	   else if (imn==0)
	   {if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]]);}
	    else         {valse=sqrt(matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1]);}
		valef=effest[imnk];
	   }
	   else if (imnk==0)
	   {if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]]);}
	    else         {valse=sqrt(matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1]);}
		valef=-effest[imn];
	   }
	  }
      else if (itp[imn]==0)
	  {if (nitp[imn]==-2)
	   {if (nitp[imnk]>0)
	    {valef=effest[imnk];
		 if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]]);}
		 else {valse=sqrt(matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1]);}
	    }
        else if (nitp[imnk]!=-1) {valef=0;valse=0;}
		else  {valef=0;valse=-1;}
	   }
	   else if ((nitp[imn]>-1) && (nitp[imnk]==-2))
	   {valef=-effest[imn];
        if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]]);}
		else {valse=sqrt(matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1]);}
	   }
	   else if ( (nitp[imn]>-1) && (nitp[imnk]==-1))
	   {valef=0;valse=-1;}
	   else if ( (nitp[imn]==-1) && (effest[imn]==0))
       {valef=0;valse=-1;
	   }
	   else if ((nitp[imnk]==-1) && (nitp[imn]==-1))
       {valef=0;valse=-1;
	   }
       else if (nitp[imn]==0)
	   {valef=effest[imnk];
        if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]]);}
	    else         {valse=sqrt(matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1]);}
	   }
	   else if (nitp[imn]>0)
	   {if ((chxt!=3) && (chxt!=6) && (chxt!=4))
        {valef= matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]];
	     valef+=matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]];
	     valef-=2*matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imnk]];
        }
		else
		{valef= matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1];
	     valef+=matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1];
	     valef-=2*matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imnk]-1];
        }
		valse=sqrt(valef);
        valef=effest[imnk]-effest[imn];
	   }
	   else
	   {valef=effest[imnk]-effest[imn];valse=0;
	   }

	  }
      else if (itp[imnk]==0)
	  {
       if (nitp[imnk]==-2)
	   {if(nitp[imn]>0)
	    {valef=-effest[imn];
	     if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]]);}
		 else {valse=sqrt(matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1]);}
        }
	    else if (nitp[imn]!=-1) {valef=0;valse=0;}
		else  {valef=0;valse=-1;}
	   }
       else if ( (nitp[imnk]>-1) && (nitp[imn]==-1))
	   {valef=0;valse=-1;}
	   else if ((nitp[imnk]==-1) && (effest[imnk]==0))
       {valef=0;valse=-1;
	   }
	   else if ((nitp[imnk]==-1) && (nitp[imn]==-1))
       {valef=0;valse=-1;
	   }
	   else if (nitp[imnk]==0)
	   {valef=-effest[imn];
        if ((chxt!=3) && (chxt!=6) && (chxt!=4)) {valse=sqrt(matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]]);}
		else  {valse=sqrt(matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1]);}
	   }
	   else if (nitp[imnk]>0)
	   {if ((chxt!=3) && (chxt!=6) && (chxt!=4))
		{valef= matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imn]];
	     valef+=matse[nbhest-1+nitp[imnk]][nbhest-1+nitp[imnk]];
	     valef-=2*matse[nbhest-1+nitp[imn]][nbhest-1+nitp[imnk]];
        }
		else
		{valef= matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imn]-1];
	     valef+=matse[nbhest-1+nitp[imnk]-1][nbhest-1+nitp[imnk]-1];
	     valef-=2*matse[nbhest-1+nitp[imn]-1][nbhest-1+nitp[imnk]-1];
        }
		valse=sqrt(valef);
        valef=effest[imnk]-effest[imn];
	   }
	   else
	   {valef=effest[imnk]-effest[imn];valse=0;
	   }
      }
      affich3(outres,valef,valse);
	  affich3(outfile,valef,valse);
      fprintf(outfile,"</td></tr>\n");fprintf(outres,"\n");
     }
     vect1=vect1->down;
    }
   }
   fprintf(outfile,"</table>\n");
  }
 }

 if ((haplozero==0) && (chxt==5))
 {
  fprintf(outfile,"<table border=0  width=80%%>\n");
  for (nc=0;nc<nkat;nc++)
  {fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
   fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
   fprintf(outfile,"<tr><td colspan=""5"">Information for CLASS LEVEL %d </td></tr>\n",nc+2);
   fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
   fprintf(outres,"\n");
   fprintf(outres,"\n");
   fprintf(outres,"Information for CLASS LEVEL %d\n",nc+2);
   fprintf(outres,"\n");

   for (i=1;i<nbloci+1;i++)
   {fprintf(outfile,"<tr><td colspan=""5"">  </td></tr>\n");
    fprintf(outfile,"<tr><td align=left colspan=""5"">Polymorphism %d %c/%c</td></tr>\n",i,letter[i-1][0],letter[i-1][1]);
    fprintf(outres,"\n");
    fprintf(outres,"Polymorphism %d %c/%c\n",i,letter[i-1][0],letter[i-1][1]);
	k=pow(2,nbloci-i);
    j=0;
    vect1=tnbhbase;trouve=0;
    while ( (vect1!=NULL) && (trouve==0))
    {imn=coding(fcoda2[vect1->numnew]); imnk=coding(fcoda2[vect1->numnew+k]);
	 if ((vect1->listall[i-1]==1) && (imn>-1) && (imnk>-1) && (freqkat[nc+1][fcoda2[vect1->numnew]]>0) &&  (freqkat[nc+1][fcoda2[vect1->numnew+k]]>0))
     {fprintf(outfile,"<tr><td colspan=""2"" align=center>Haplotypic Background</td>\n");
      fprintf(outfile,"<td align=left>");
	  fprintf(outres,"Haplotypic Background ");
	  for (ii=0;ii<nbloci;ii++)
      {if (ii+1!=i)
       {l=letter[ii][0]*(vect1->listall[ii]==1)+letter[ii][1]*(vect1->listall[ii]==2);
        fprintf(outfile,"%c",l);fprintf(outres,"%c",l);
	   }
       else if (ii+1==i) {fprintf(outfile,"-");
	                      fprintf(outres,"-");
					      }
      }
	  fprintf(outfile,"</td><td align= left width=10%%>");
	  fprintf(outfile," ID %d-%d \t",imn*nkat+nc,imnk*nkat+nc);
	  fprintf(outfile,"</td><td colspan=""2"">");
	  fprintf(outres,"\tID %d-%d\t",imn*nkat+nc,imnk*nkat+nc);

	  idn=imn*nkat+nc;
	  idk=imnk*nkat+nc;

      if ((idn==nc) && (itp[idk]==1))
	  {valef=effest[idk];
	   valse=sqrt(matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]]);
	  }
      else if ((idn==nc) && (itp[idk]==0))
	  {valef=effest[idk];
       if (valef!=0)
	   {valse=sqrt(matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]]);
	   }
       else valse=-1;
	  }
      else if ((idk==nc) && (itp[idn]==1))
	  {valef=-effest[idn];
       valse=sqrt(matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]]);
	  }
	  else  if ((itp[idn]==1) && (itp[idk]==1))
	  {if ((idn!=0) && (idk!=0))
       {valef=matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]];
	    valef+=matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]];
	    valef-=2*matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idk]];
        valse=sqrt(valef);
        valef=effest[idk]-effest[idn];
	   }
	   else if (idn==0)
	   {valse=sqrt(matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]]);
	   	valef=effest[idk];
	   }
	   else if (idk==0)
	   {valse=sqrt(matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]]);
	    valef=-effest[idn];
	   }
	  }
	  else if (itp[idn]==0)
	  {if (nitp[idn]==-2)
	   {if (nitp[idk]>0)
	    {valef=effest[idk];valse=sqrt(matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]]);
		}
        else if (nitp[idk]!=-1) {valef=0;valse=0;}
		else  {valef=0;valse=-1;}
	   }
	   else if ((nitp[idn]>-1) && (nitp[idk]==-2))
	   {valef=-effest[idn];valse=sqrt(matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]]);
	   }
	   else if ( (nitp[idn]>-1) && (nitp[idk]==-1))
	   {valef=0;valse=-1;}
	   else if ( (nitp[idn]==-1) && (effest[idn]==0))
       {valef=0;valse=-1;}
	   else if ((nitp[idk]==-1) && (nitp[idn]==-1))
       {valef=0;valse=-1;}
       else if (nitp[idn]==0)
	   {valef=effest[idk];valse=sqrt(matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]]);}
	   else if (nitp[idn]>0)
	   {valef= matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]];
	    valef+=matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]];
	    valef-=2*matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idk]];
       	valse=sqrt(valef);
        valef=effest[idk]-effest[idn];
	   }
	   else
	   {valef=effest[idk]-effest[idn];valse=0;
	   }
	  }
      else if (itp[idk]==0)
	  {if (nitp[idk]==-2)
	   {if(nitp[idn]>0)
	    {valef=-effest[idn];valse=sqrt(matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]]);}
	    else if (nitp[idn]!=-1) {valef=0;valse=0;}
		else  {valef=0;valse=-1;}
	   }
       else if ( (nitp[idk]>-1) && (nitp[idn]==-1))
	   {valef=0;valse=-1;}
	   else if ((nitp[idk]==-1) && (effest[idk]==0))
       {valef=0;valse=-1;
	   }
	   else if ((nitp[idk]==-1) && (nitp[idn]==-1))
       {valef=0;valse=-1;
	   }
	   else if (nitp[idk]==0)
	   {valef=-effest[idn];valse=sqrt(matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]]);}
	   /*else if (nitp[idk]>0)
	   {valef= matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idn]];
	    valef+=matse[nbhest-1+nitp[idk]][nbhest-1+nitp[idk]];
	    valef-=2*matse[nbhest-1+nitp[idn]][nbhest-1+nitp[idk]];
        valse=sqrt(valef);
        valef=effest[idk]-effest[idn];
	   }*/
	   else
	   {valef=effest[idk]-effest[idn];valse=0;}
      }
	  else {valef=0;valse=-1;}


	  affich3(outfile,valef,valse);
	  affich3(outres,valef,valse);
	  fprintf(outfile,"</td></tr>\n");fprintf(outres,"\n");

     }
	 vect1=vect1->down;
    }
   }/*i*/
  }/*nc*/
  fprintf(outfile,"</table>\n");
 }
//printf("BOU\n");
 if (chxt==2)  {phenomean(outres,outfile,matse);}
 //printf("BOU3\n");


 if (rsq==1)
 {fprintf(outres,"\n\n\n");
  fprintf(outfile,"<table border=0  width=80%%>\n");
  if (xlnk==0) rsquare(moyfreq,outfile,outres);
  fprintf(outfile,"</table>");
 }

 if ((ldmatrix==1) && (xlnk==1))
 {fprintf(outres,"\n\n\n");
  fprintf(outfile,"<table border=0  width=80%%>\n");
  ldp=(double **) malloc ((size_t) (nbloci*sizeof(double *)));
  lr2=(double **) malloc ((size_t) (nbloci*sizeof(double *)));
  chidp=(double **) malloc ((size_t) (nbloci*sizeof(double *)));
  for (i=0;i<nbloci;i++)
  {ldp[i]=(double *) malloc((size_t) (nbloci*sizeof(double)));
   lr2[i]=(double *) malloc((size_t) (nbloci*sizeof(double)));
   chidp[i]=(double *) malloc((size_t) (nbloci*sizeof(double)));
  }
  for (k=0;k<nbloci;k++) for (l=0;l<nbloci;l++)
  {chidp[k][l]=0;
   if (k==l) {ldp[k][l]=1;lr2[k][l]=1;}
   else {ldp[k][l]=0;lr2[k][l]=0;}
  }

  if ((chxt==0) || (chxt==2))
  {for (k=0;k<nbloci;k++) alfreq[k]=0;
   vect1=tnbhbase;
   while (vect1!=NULL)
   {ll=fcoda2[vect1->numnew];
    if ((vect1->present==1) && (moyfreq[ll]>0))
    {for (k=0;k<nbloci;k++) alfreq[k]+=moyfreq[ll]*(vect1->listall[k]==1);
    }
    vect1=vect1->down;
   }
   vect1=tnbhbase;
   while (vect1!=NULL)
   {ll=fcoda2[vect1->numnew];
    if ((vect1->present==1) && (moyfreq[ll]>0))
    {for (k=0;k<nbloci;k++)
     {ik=2;if (alfreq[k]<0.5) {ik=1;}
      for (ii=k+1;ii<nbloci;ii++)
      {il=2;if (alfreq[ii]<0.5) {il=1;}
       if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
       {ldp[k][ii]+=moyfreq[ll];
       }
      }
	 }
    }
    vect1=vect1->down;
   }
   for (k=0;k<nbloci;k++)
   for (ii=k+1;ii<nbloci;ii++)
   {mink=alfreq[k]*(alfreq[k]<0.5)+(1-alfreq[k])*(alfreq[k]>=0.5);
    minl=alfreq[ii]*(alfreq[ii]<0.5)+(1-alfreq[ii])*(alfreq[ii]>=0.5);
	ldp[k][ii]-=mink*minl;
	chidp[k][ii]=(ldp[k][ii]*ldp[k][ii]*(2*nbhf[1][0]+nbhf[0][0]))/(alfreq[k]*alfreq[ii]*(1-alfreq[k])*(1-alfreq[ii]));
	if (ldp[k][ii]<0) {dmax=mink*minl;}
	else
	{if (mink<minl)	{dmax=mink*(1-minl);}
     else {dmax=minl*(1-mink);}
    }
    lr2[k][ii]=ldp[k][ii]*ldp[k][ii]/(mink*minl*(1-mink)*(1-minl));lr2[ii][k]=lr2[k][ii];
    ldp[k][ii]/=dmax;
    if (fabs(ldp[k][ii])>1.0) {ldp[k][ii]/=fabs(ldp[k][ii]);}
   }

   fprintf(outfile,"<tr><td width=60%%> </td><td width=40%%> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX</td></tr>\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outres,"\t\t\tLD MATRIX\n\n");

   for (k=0;k<nbloci;k++)
   {fprintf(outfile,"<tr><td align=right>");
    fprintf(outres,"\t\t");
    for (ii=0;ii<nbloci;ii++)
    {if (ii==k) {fprintf(outfile," 1");fprintf(outres," 1\t");}
     else if (ii<k) {fprintf(outfile," %3.2f",lr2[k][ii]);
                     fprintf(outres," %3.2f\t",lr2[k][ii]);
                    }
     else if (ii>k) {fprintf(outfile," %3.2f",ldp[k][ii]);
                     fprintf(outres," %3.2f\t",ldp[k][ii]);
                    }
    }
    fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
   }


   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outres,"\n\n");
   fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value MATRIX</td></tr>\n");
   fprintf(outres,"\t\t\tP-Value MATRIX\n\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   for (k=0;k<nbloci;k++)
   {fprintf(outfile,"<tr><td align=right>");fprintf(outres,"\t\t");
     fprintf(outfile," X ");fprintf(outres," X\t");
    for (ii=k+1;ii<nbloci;ii++) {fprintf(outfile," %1.3f",chdtrc(1.0,chidp[k][ii]));
	                             fprintf(outres," %1.3f\t",chdtrc(1.0,chidp[k][ii]));
								}
    fprintf(outfile,"</td><td> </td></tr>");
	fprintf(outres,"\n");for (ii=0;ii<k+1;ii++) fprintf(outres,"\t");
   }
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outfile,"\n\n");
  }
  else if (chxt==1)
  {for (i=0;i<2;i++)
   {for (k=0;k<nbloci;k++)
     for (ii=0;ii<nbloci;ii++)
      {chidp[k][ii]=0;
       if (k==ii) {ldp[k][ii]=1;lr2[k][ii]=1;}
        else {ldp[k][ii]=0;lr2[k][ii]=0;}
      }
    for (k=0;k<nbloci;k++) alfreq[k]=0;
    vect1=tnbhbase;
    if (i==0)
	{while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqtem[ll]>0))
      {for (k=0;k<nbloci;k++) alfreq[k]+=freqtem[ll]*(vect1->listall[k]==1);
      }
      vect1=vect1->down;
     }
     vect1=tnbhbase;
     while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqtem[ll]>0))
      {for (k=0;k<nbloci;k++)
       {ik=2;if (alfreq[k]<0.5) {ik=1;}
        for (ii=k+1;ii<nbloci;ii++)
        {il=2;if (alfreq[ii]<0.5) {il=1;}
         if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
         {ldp[k][ii]+=freqtem[ll];
         }
        }
	   }
      }
      vect1=vect1->down;
     }
   	}
    if (i==1)
	{while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqcas[ll]>0))
      {for (k=0;k<nbloci;k++) alfreq[k]+=freqcas[ll]*(vect1->listall[k]==1);
      }
      vect1=vect1->down;
     }
     vect1=tnbhbase;
     while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqcas[ll]>0))
      {for (k=0;k<nbloci;k++)
       {ik=2;if (alfreq[k]<0.5) {ik=1;}
        for (ii=k+1;ii<nbloci;ii++)
        {il=2;if (alfreq[ii]<0.5) {il=1;}
         if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
         {ldp[k][ii]+=freqcas[ll];
         }
        }
	   }
      }
      vect1=vect1->down;
     }
   	}
    for (k=0;k<nbloci;k++)
    for (ii=k+1;ii<nbloci;ii++)
    {mink=alfreq[k]*(alfreq[k]<0.5)+(1-alfreq[k])*(alfreq[k]>=0.5);
     minl=alfreq[ii]*(alfreq[ii]<0.5)+(1-alfreq[ii])*(alfreq[ii]>=0.5);
     ldp[k][ii]-=mink*minl;
     chidp[k][ii]=(ldp[k][ii]*ldp[k][ii])/(alfreq[k]*alfreq[ii]*(1-alfreq[k])*(1-alfreq[ii]));
     if (i==0) {chidp[k][ii]*=(2*nbhf[1][1]+nbhf[0][1]);}
     else if (i==1) {chidp[k][ii]*=(2*nbhf[1][2]+nbhf[0][2]);}
	 if (ldp[k][ii]<0) {dmax=mink*minl;}
	 else
	 {if (mink<minl)	{dmax=mink*(1-minl);}
      else {dmax=minl*(1-mink);}
     }
     lr2[k][ii]=ldp[k][ii]*ldp[k][ii]/(mink*minl*(1-mink)*(1-minl));lr2[ii][k]=lr2[k][ii];
     ldp[k][ii]/=dmax;
     if (fabs(ldp[k][ii])>1.0) {ldp[k][ii]/=fabs(ldp[k][ii]);}
    }

    fprintf(outres,"\n\n");
    fprintf(outfile,"<tr><td width=60%%> </td></tr><tr><td width=40%%> </td></tr>\n");
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    if (i==0)
    {fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX IN CONTROLS</td></tr>\n");
     fprintf(outres,"\t\t\tLD MATRIX IN CONTROLS\n\n");
    }
    if (i==1)
    {fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX IN CASES</td></tr>\n");
     fprintf(outres,"\t\t\tLD MATRIX IN CASES\n\n");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    for (k=0;k<nbloci;k++)
    {fprintf(outres,"\t\t");
     fprintf(outfile,"<tr><td align=right>");
     for (ii=0;ii<nbloci;ii++)
     {if (ii==k) {fprintf(outfile," 1");fprintf(outres," 1\t");}
      else if (ii<k) {fprintf(outfile," %3.2f",lr2[k][ii]);
                      fprintf(outres," %3.2f\t",lr2[k][ii]);
                     }
      else if (ii>k) {fprintf(outfile," %3.2f",ldp[k][ii]);
                      fprintf(outres," %3.2f\t",ldp[k][ii]);
                     }
     }
     fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n"); fprintf(outres,"\n\n");
    if (i==0)
	{fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value IN CONTROLS</td></tr>\n");
     fprintf(outres,"\t\t\tP-Value IN CONTROLS\n\n");
	}
    if (i==1)
	{fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value IN CASES</td></tr>\n");
     fprintf(outres,"\t\t\tP-Value IN CASES\n\n");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   	for (k=0;k<nbloci;k++)
    {fprintf(outfile,"<tr><td align=right>");fprintf(outres,"\t\t");
	 fprintf(outfile," X ");fprintf(outres," X\t"); 
       for (ii=k+1;ii<nbloci;ii++) {fprintf(outfile,"%.4f\t",chdtrc(1.0,chidp[k][ii]));
	                              fprintf(outres,"%.4f\t",chdtrc(1.0,chidp[k][ii]));
								 }
     fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
	 for (ii=0;ii<k+1;ii++) fprintf(outres,"\t");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outres,"\n");
   }
  }
  fprintf(outfile,"</table></body></html>");
 }
 else if ((ldmatrix==1) && (xlnk==0))
 {fprintf(outres,"\n\n\n");
  fprintf(outfile,"<table border=0  width=80%%>\n");


  ldp=(double **) malloc ((size_t) (nbloci*sizeof(double *)));
  lr2=(double **) malloc ((size_t) (nbloci*sizeof(double *)));
  chidp=(double **) malloc ((size_t) (nbloci*sizeof(double *)));
  for (i=0;i<nbloci;i++)
  {ldp[i]=(double *) malloc((size_t) (nbloci*sizeof(double)));
   lr2[i]=(double *) malloc((size_t) (nbloci*sizeof(double)));
   chidp[i]=(double *) malloc((size_t) (nbloci*sizeof(double)));
  }
  for (k=0;k<nbloci;k++)
   for (l=0;l<nbloci;l++)
	{chidp[k][l]=0;
     if (k==l) {ldp[k][l]=1;lr2[k][l]=1;} else {ldp[k][l]=0;lr2[k][l]=0;}
	}

  if ((chxt==0) || (chxt==2))
  {for (k=0;k<nbloci;k++) alfreq[k]=0;
   vect1=tnbhbase;
   while (vect1!=NULL)
   {ll=fcoda2[vect1->numnew];
    if ((vect1->present==1) && (moyfreq[ll]>0))
    {for (k=0;k<nbloci;k++) alfreq[k]+=moyfreq[ll]*(vect1->listall[k]==1);
    }
    vect1=vect1->down;
   }
   vect1=tnbhbase;
   while (vect1!=NULL)
   {ll=fcoda2[vect1->numnew];
    if ((vect1->present==1) && (moyfreq[ll]>0))
    {for (k=0;k<nbloci;k++)
     {ik=2;if (alfreq[k]<0.5) {ik=1;}
      for (ii=k+1;ii<nbloci;ii++)
      {il=2;if (alfreq[ii]<0.5) {il=1;}
       if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
       {ldp[k][ii]+=moyfreq[ll];
       }
      }
	 }
    }
    vect1=vect1->down;
   }
   for (k=0;k<nbloci;k++)
   for (ii=k+1;ii<nbloci;ii++)
   {mink=alfreq[k]*(alfreq[k]<0.5)+(1-alfreq[k])*(alfreq[k]>=0.5);
    minl=alfreq[ii]*(alfreq[ii]<0.5)+(1-alfreq[ii])*(alfreq[ii]>=0.5);
	ldp[k][ii]-=mink*minl;
	chidp[k][ii]=2*(ldp[k][ii]*ldp[k][ii]*nbused)/(alfreq[k]*alfreq[ii]*(1-alfreq[k])*(1-alfreq[ii]));
	if (ldp[k][ii]<0) {dmax=mink*minl;}
	else
	{if (mink<minl)	{dmax=mink*(1-minl);}
     else {dmax=minl*(1-mink);}
    }
    lr2[k][ii]=ldp[k][ii]*ldp[k][ii]/(mink*minl*(1-mink)*(1-minl));lr2[ii][k]=lr2[k][ii];
    ldp[k][ii]/=dmax;
    if (fabs(ldp[k][ii])>1.0) {ldp[k][ii]/=fabs(ldp[k][ii]);}
   }

   fprintf(outfile,"<tr><td width=60%%> </td><td width=40%%> </td></tr>\n");
   fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX</td></tr>\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   fprintf(outres,"\t\t\tLD MATRIX\n\n");

   for (k=0;k<nbloci;k++)
   {fprintf(outfile,"<tr><td align=right>");fprintf(outres,"\t\t");
    for (ii=0;ii<nbloci;ii++)
    {if (ii==k) {fprintf(outfile," 1");fprintf(outres," 1\t");}
     else if (ii<k) {fprintf(outfile," %3.2f",lr2[k][ii]);
                     fprintf(outres," %3.2f\t",lr2[k][ii]);
                    }
     else if (ii>k) {fprintf(outfile," %3.2f",ldp[k][ii]);
                     fprintf(outres," %3.2f\t",ldp[k][ii]);
                    }
    }
    fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
   }
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outres,"\n\n");
   fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value MATRIX</td></tr>\n");
   fprintf(outres,"\t\t\tP-Value MATRIX\n\n");
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   for (k=0;k<nbloci;k++)
   {fprintf(outfile,"<tr><td align=right>");fprintf(outres,"\t\t");
     fprintf(outfile," X ");fprintf(outres," X\t");
     for (ii=k+1;ii<nbloci;ii++) {fprintf(outfile," %1.3f",chdtrc(1.0,chidp[k][ii]));
	                             fprintf(outres," %1.3f\t",chdtrc(1.0,chidp[k][ii]));
								}
    fprintf(outfile,"</td><td> </td></tr>");
	fprintf(outres,"\n");for (ii=0;ii<k+1;ii++) fprintf(outres,"\t");
   }
   fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outfile,"\n\n");
  }
  else if ((chxt==1) || (chxt==3))
  {for (i=0;i<2;i++)
   {for (k=0;k<nbloci;k++)
     for (ii=0;ii<nbloci;ii++)
	  {chidp[k][ii]=0;
       if (k==ii) {ldp[k][ii]=1;lr2[k][ii]=1;}
	   else {ldp[k][ii]=0;lr2[k][ii]=0;}
      }
    for (k=0;k<nbloci;k++) alfreq[k]=0;
    vect1=tnbhbase;
    if (i==0)
	{while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqtem[ll]>0))
      {for (k=0;k<nbloci;k++) alfreq[k]+=freqtem[ll]*(vect1->listall[k]==1);
      }
      vect1=vect1->down;
     }
     vect1=tnbhbase;
     while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqtem[ll]>0))
      {for (k=0;k<nbloci;k++)
       {ik=2;if (alfreq[k]<0.5) {ik=1;}
        for (ii=k+1;ii<nbloci;ii++)
        {il=2;if (alfreq[ii]<0.5) {il=1;}
         if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
         {ldp[k][ii]+=freqtem[ll];
         }
        }
	   }
      }
      vect1=vect1->down;
     }
   	}
    if (i==1)
	{while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqcas[ll]>0))
      {for (k=0;k<nbloci;k++) alfreq[k]+=freqcas[ll]*(vect1->listall[k]==1);
      }
      vect1=vect1->down;
     }
     vect1=tnbhbase;
     while (vect1!=NULL)
     {ll=fcoda2[vect1->numnew];
      if ((vect1->present==1) && (freqcas[ll]>0))
      {for (k=0;k<nbloci;k++)
       {ik=2;if (alfreq[k]<0.5) {ik=1;}
        for (ii=k+1;ii<nbloci;ii++)
        {il=2;if (alfreq[ii]<0.5) {il=1;}
         if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
         {ldp[k][ii]+=freqcas[ll];
         }
        }
	   }
      }
      vect1=vect1->down;
     }
   	}
    for (k=0;k<nbloci;k++)
    for (ii=k+1;ii<nbloci;ii++)
    {mink=alfreq[k]*(alfreq[k]<0.5)+(1-alfreq[k])*(alfreq[k]>=0.5);
     minl=alfreq[ii]*(alfreq[ii]<0.5)+(1-alfreq[ii])*(alfreq[ii]>=0.5);
     ldp[k][ii]-=mink*minl;
     chidp[k][ii]=2*(ldp[k][ii]*ldp[k][ii])/(alfreq[k]*alfreq[ii]*(1-alfreq[k])*(1-alfreq[ii]));
     if (i==0) {chidp[k][ii]*=nbtem;}
     else if (i==1) {chidp[k][ii]*=nbcas;}
	 if (ldp[k][ii]<0) {dmax=mink*minl;}
	 else
	 {if (mink<minl)	{dmax=mink*(1-minl);}
      else {dmax=minl*(1-mink);}
     }
     lr2[k][ii]=ldp[k][ii]*ldp[k][ii]/(mink*minl*(1-mink)*(1-minl));lr2[ii][k]=lr2[k][ii];
     ldp[k][ii]/=dmax;
     if (fabs(ldp[k][ii])>1.0) {ldp[k][ii]/=fabs(ldp[k][ii]);}
    }

    fprintf(outres,"\n\n");
	fprintf(outfile,"<tr><td width=60%%> </td></tr><tr><td width=40%%> </td></tr>\n");
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    if (i==0)
	{fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX IN CONTROLS</td></tr>\n");
     fprintf(outres,"\t\t\tLD MATRIX IN CONTROLS\n\n");
    }
    if (i==1)
	{fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX IN CASES</td></tr>\n");
     fprintf(outres,"\t\t\tLD MATRIX IN CASES\n\n");
    }
     fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    for (k=0;k<nbloci;k++)
    {fprintf(outres,"\t\t");
     fprintf(outfile,"<tr><td align=right>");

     for (ii=0;ii<nbloci;ii++)
     {if (ii==k) {fprintf(outfile," 1");fprintf(outres," 1\t");}
      else if (ii<k) {fprintf(outfile," %3.2f",lr2[k][ii]);
                      fprintf(outres," %3.2f\t",lr2[k][ii]);
                     }
      else if (ii>k) {fprintf(outfile," %3.2f",ldp[k][ii]);
                      fprintf(outres," %3.2f\t",ldp[k][ii]);
                     }
     }
    fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n"); fprintf(outres,"\n\n");
    if (i==0)
	{fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value IN CONTROLS</td></tr>\n");
     fprintf(outres,"\t\t\tP-Value IN CONTROLS\n\n");
	}
    if (i==1)
	{fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value IN CASES</td></tr>\n");
     fprintf(outres,"\t\t\tP-Value IN CASES\n\n");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   	for (k=0;k<nbloci;k++)
    {fprintf(outfile,"<tr><td align=right>");fprintf(outres,"\t\t");
     fprintf(outfile," X ");fprintf(outres," X\t");	
     for (ii=k+1;ii<nbloci;ii++) {fprintf(outfile,"%.4f\t",chdtrc(1.0,chidp[k][ii]));
	                              fprintf(outres,"%.4f\t",chdtrc(1.0,chidp[k][ii]));
								 }
     fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
	 for (ii=0;ii<k+1;ii++) fprintf(outres,"\t");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outres,"\n");
   }
  }
  else if (chxt==5)
  {for (i=0;i<nbcatego;i++)
   {for (k=0;k<nbloci;k++)
     for (ii=0;ii<nbloci;ii++)
	  {chidp[k][ii]=0;
       if (k==ii) {ldp[k][ii]=1;lr2[k][ii]=1;}
	   else {ldp[k][ii]=0;lr2[k][ii]=0;}
      }
    for (k=0;k<nbloci;k++) alfreq[k]=0;
    vect1=tnbhbase;

	while (vect1!=NULL)
    {ll=fcoda2[vect1->numnew];
     if ((vect1->present==1) && (freqkat[i][ll]>0))
     {for (k=0;k<nbloci;k++) alfreq[k]+=freqkat[i][ll]*(vect1->listall[k]==1);
     }
     vect1=vect1->down;
    }
    vect1=tnbhbase;
    while (vect1!=NULL)
    {ll=fcoda2[vect1->numnew];
     if ((vect1->present==1) && (freqkat[i][ll]>0))
     {for (k=0;k<nbloci;k++)
      {ik=2;if (alfreq[k]<0.5) {ik=1;}
       for (ii=k+1;ii<nbloci;ii++)
       {il=2;if (alfreq[ii]<0.5) {il=1;}
        if ((vect1->listall[k]==ik) && (vect1->listall[ii]==il))
        {ldp[k][ii]+=freqkat[i][ll];
        }
       }
	  }
     }
     vect1=vect1->down;
    }
   	for (k=0;k<nbloci;k++)
    for (ii=k+1;ii<nbloci;ii++)
    {mink=alfreq[k]*(alfreq[k]<0.5)+(1-alfreq[k])*(alfreq[k]>=0.5);
     minl=alfreq[ii]*(alfreq[ii]<0.5)+(1-alfreq[ii])*(alfreq[ii]>=0.5);
     ldp[k][ii]-=mink*minl;
     chidp[k][ii]=2*(ldp[k][ii]*ldp[k][ii])/(alfreq[k]*alfreq[ii]*(1-alfreq[k])*(1-alfreq[ii]));
     chidp[k][ii]*=nbsujktgo[i];
	 if (ldp[k][ii]<0) {dmax=mink*minl;}
	 else
	 {if (mink<minl)	{dmax=mink*(1-minl);}
      else {dmax=minl*(1-mink);}
     }
     lr2[k][ii]=ldp[k][ii]*ldp[k][ii]/(mink*minl*(1-mink)*(1-minl));lr2[ii][k]=lr2[k][ii];
     ldp[k][ii]/=dmax;
     if (fabs(ldp[k][ii])>1.0) {ldp[k][ii]/=fabs(ldp[k][ii]);}
    }

    fprintf(outres,"\n\n");
	fprintf(outfile,"<tr><td width=60%%> </td></tr><tr><td width=40%%> </td></tr>\n");
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
    fprintf(outfile,"<tr><td align=center colspan=""2"">LD MATRIX IN CLASS LEVEL %d</td></tr>\n",i+1);
    fprintf(outres,"\t\t\tLD MATRIX IN CLASS LEVEL %d\n\n",i+1);
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");

    for (k=0;k<nbloci;k++)
    {fprintf(outres,"\t\t");
     fprintf(outfile,"<tr><td align=right>");
     for (ii=0;ii<nbloci;ii++)
     {if (ii==k) {fprintf(outfile," 1");fprintf(outres," 1\t");}
      else if (ii<k) {fprintf(outfile," %3.2f",lr2[k][ii]);
                      fprintf(outres," %3.2f\t",lr2[k][ii]);
                     }
      else if (ii>k) {fprintf(outfile," %3.2f",ldp[k][ii]);
                      fprintf(outres," %3.2f\t",ldp[k][ii]);
                     }
     }
    fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n"); fprintf(outres,"\n\n");

    fprintf(outfile,"<tr><td align=center colspan=""2"">P-Value IN CLASS LEVEL %d</td></tr>\n",i+1);
    fprintf(outres,"\t\t\tP-Value IN CLASS LEVEL %d\n\n",i+1);

    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");
   	for (k=0;k<nbloci;k++)
    {fprintf(outfile,"<tr><td align=right>");fprintf(outres,"\t\t");
     fprintf(outfile," X ");fprintf(outres," X\t");
     for (ii=k+1;ii<nbloci;ii++) {fprintf(outfile,"%.4f\t",chdtrc(1.0,chidp[k][ii]));
	                              fprintf(outres,"%.4f\t",chdtrc(1.0,chidp[k][ii]));
								 }
     fprintf(outfile,"</td><td> </td></tr>");fprintf(outres,"\n");
	 for (ii=0;ii<k+1;ii++) fprintf(outres,"\t");
    }
    fprintf(outfile,"<tr><td> </td></tr><tr><td> </td></tr>\n");fprintf(outres,"\n");
   }
  }
  fprintf(outfile,"</table>");
 }

/* POUR JAVA CECI EST INDISPONIBLE
 if ((chxt>0) && (haplozero==0) && (hypoth==0) && (interor==0) && (chxt<5))
 {wald=0;
  printf("Do you want some variance-covariance matrix to be display (y/n)? \n");
  do
  {scanf("%c",&rep);
  }while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));
  if ((rep=='y')  || (rep=='Y')) {wald=1;}

  if (wald==1)
  {printf("How many parameters involved in this matrix ?\n");
   do
   {scanf("%d",&nbwald);
   }while (nbwald<2);

   idwald =(int *)malloc((size_t) (nbwald*sizeof(int)));

   for (i=0;i<nbwald;i++)
   {do
    {printf("Please enter the number of haplotypes %d : ",i+1);
     scanf("%d",&idwald[i]);
     if (itp[idwald[i]-1]==0) {printf("Choice impossible because this effect was not estimated\n");}
    }while ( (idwald[i]<2) || (itp[idwald[i]-1]==0));
   }

   fprintf(outres,"\n\nInformation for Generalized Test Statistic\nVector of regression parameters:\n");
   for (i=0;i<nbwald;i++) fprintf(outres,"Haplotype effect %d : %lf\n",idwald[i],effest[idwald[i]-1]);

   fprintf(outres,"\n\nVariance-Covariance matrix of regression parameters:\n");
   for (i=0;i<nbwald;i++)
   {for (j=0;j<nbwald;j++)
	fprintf(outres,"%lf  ",matse[nbhest-1+nitp[idwald[i]-1]][nbhest-1+nitp[idwald[j]-1]]);
	fprintf(outres,"\n");
   }


   fprintf(outfile,"</br><table>");
   fprintf(outfile,"<tr><td>Information for Generalized Test Statistic\nVector of regression parameters</td></tr>");
   for (i=0;i<nbwald;i++)
	 fprintf(outfile,"<tr><td>Haplotype effect %d :&nbsp; %lf</td></tr>\n",idwald[i],effest[idwald[i]-1]);

   fprintf(outfile,"<tr></tr><tr></tr>");
   fprintf(outfile,"<tr><td>Variance-Covariance matrix of regression parameters </td></tr>\n");

   for (i=0;i<nbwald;i++)
   {fprintf(outfile,"<tr><td>");
    for (j=0;j<nbwald;j++) fprintf(outfile,"%.6f  &nbsp;&nbsp; ",matse[nbhest-1+nitp[idwald[i]-1]][nbhest-1+nitp[idwald[j]-1]]);
    fprintf(outfile,"</td></tr>");
   }
   fprintf(outfile,"</table>");

   idwald=NULL;
   free(idwald);
  }
 }
 if ((chxt==5) && (haplozero==0) && (hypoth==0) && (interor==0))
 {wald=0;
  printf("Do you want some variance-covariance matrix to be display (y/n)? \n");
  do
  {scanf("%c",&rep);
  }while ( (rep!='y') && (rep!='Y') && (rep!='n') && (rep!='N'));
  if ((rep=='y')  || (rep=='Y')) {wald=1;}

  if (wald==1)
  {printf("How many parameters involved in this matrix ?\n");
   do
   {scanf("%d",&nbwald);
   }while (nbwald<2);

   idwald =(int *)malloc((size_t) (nbwald*sizeof(int)));

   for (i=0;i<nbwald;i++)
   {do
    {printf("Please enter the ID number of haplotypes %d : ",i+1);
     scanf("%d",&idwald[i]);
     if (itp[idwald[i]]==0) {printf("Choice impossible because this effect was not estimated\n");}
    }while ( (idwald[i]<0) || (itp[idwald[i]]==0));
   }

   fprintf(outres,"\n\nInformation for Generalized Test Statistic\nVector of regression parameters:\n");
   for (i=0;i<nbwald;i++) fprintf(outres,"Haplotype effect %d : %lf\n",idwald[i],effest[idwald[i]]);

   fprintf(outres,"\n\nVariance-Covariance matrix of regression parameters:\n");
   for (i=0;i<nbwald;i++)
   {for (j=0;j<nbwald;j++)
	fprintf(outres,"%lf  ",matse[nbhest-1+nitp[idwald[i]]][nbhest-1+nitp[idwald[j]]]);
	fprintf(outres,"\n");
   }


   fprintf(outfile,"</br><table>");
   fprintf(outfile,"<tr><td>Information for Generalized Test Statistic\nVector of regression parameters</td></tr>");
   for (i=0;i<nbwald;i++)
	 fprintf(outfile,"<tr><td>Haplotype effect %d :&nbsp; %lf</td></tr>\n",idwald[i],effest[idwald[i]]);

   fprintf(outfile,"<tr></tr><tr></tr>");
   fprintf(outfile,"<tr><td>Variance-Covariance matrix of regression parameters </td></tr>\n");

   for (i=0;i<nbwald;i++)
   {fprintf(outfile,"<tr><td>");
    for (j=0;j<nbwald;j++) fprintf(outfile,"%.6f  &nbsp;&nbsp; ",matse[nbhest-1+nitp[idwald[i]]][nbhest-1+nitp[idwald[j]]]);
    fprintf(outfile,"</td></tr>");
   }
   fprintf(outfile,"</table>");

   idwald=NULL;
   free(idwald);
  }
 }
*/




 fprintf(outfile,"</table></body></html>");

 // New DAvid
 if (ldmatrix==1)
 {for (i=0;i<nbloci;i++)  {free(ldp[i]);ldp[i]=NULL;free(chidp[i]);chidp[i]=NULL;free(lr2[i]);lr2[i]=NULL;}
  ldp=NULL;free (ldp); chidp=NULL; free (chidp);lr2=NULL;free (lr2);
 }

 if (outfile  != NULL) fclose( outfile);outfile=NULL;
 if (outres   != NULL) fclose( outres);outres=NULL;

 free (dmat2);dmat2=NULL;
 free (mdvs2);mdvs2=NULL;
 free (modif2);modif2=NULL;
 free(tabres);tabres=NULL;
 free(moyfreq); moyfreq=NULL;
 free (fcoda1);fcoda1=NULL;
 free (fcoda2);fcoda2=NULL;
 free (tabmq);tabmq=NULL;
  free (moyeff);moyeff=NULL;
 free (effest);effest=NULL;
 free (itp);itp=NULL;
 free (itptp);itptp=NULL;
 free (nitp);nitp=NULL;
 free (nitptp);nitptp=NULL;
  free (alfreq);alfreq=NULL;
 free (numhap);numhap=NULL;
 free (inclus);inclus=NULL;
 free (hafreq);hafreq=NULL;
  if (chxt==0)
  {free (frqobs);frqobs=NULL;
    free (frqord);frqord=NULL;
    free(numord);numord=NULL;
    free(numobs);numobs=NULL;
  }
  free (place);place=NULL;
  free (hwc);hwc=NULL;
  free (freqest); freqest=NULL;
  free (tempfreq);tempfreq=NULL;
 if ( (chxt==1) ||  (chxt==3) || (chxt==4) || (chxt==6))
 {free(freqcas); freqcas=NULL; free(freqtem); freqtem=NULL;
 }
 if (chxt==5)
 {free(freqkat);freqkat=NULL;
 }
  if (chxt==3)
  {free ((double **) tablo);tablo=NULL;
   free (tabpi); tabpi=NULL;
  } 
  if ( (chxt==3) || (chxt==4)) {free(bresl);bresl=NULL;} 
   // NEW DAVID

   free(base);base=NULL;
 free(suiv);suiv=NULL;
 free(tnbhnew);tnbhnew=NULL;
 free(tnbhbase);tnbhbase=NULL;
 free (vect1); vect1=NULL;
 free (vect2); vect2=NULL;
 

     
  /////
 //printf("DEBUT\n"); 
  if (chxt>2)
 {for (i=0;i<nall;i++)
  {free(vecbeta[i]); vecbeta[i]=NULL;free(vecbetat[i]); vecbetat[i]=NULL;
   free(vintra[i]); vintra[i]=NULL;free(vinter[i]); vinter[i]=NULL;
  }
  free(vinter);free(vintra);vinter=NULL; vintra=NULL;
  free(vecbeta);free(vecbetat);vecbeta=NULL; vecbetat=NULL;
 }
//printf("F1\n");
//printf("nbadd= %d\n",nbadd);
 if (nbadd>0)
 {for (i=0;i<nbadd;i++) {free(tadd[i]);tadd[i]=NULL;}
  free (tadd);tadd=NULL;
 }
 //printf("F2\n");
 //printf("intercov= %d\n",intercov);
 if (intercov>0)
 {for (i=0;i<intercov;i++) {free(tabint[i]);tabint[i]=NULL;}
  free (tabint);tabint=NULL;
 }
//printf("F3\n");
 //printf("chxt= %d\n",chxt);
 if ( (chxt!=0) && (chxt!=2) )
 {for (i=0;i<nbhhypo;i++)
  {free(freqdist[i]);free(tempdist[i]);freqdist[i]=NULL; tempdist[i]=NULL;}
  free (freqdist);freqdist=NULL;
  free (tempdist);tempdist=NULL;
 }
//printf("F4= %d\n",nbhypor);
//printf("INTEROR= %d\n",interor);
//printf("nbhypor= %d\n",nbhypor);
 if (nbhypor>0)
 {for (i=0;i<nbhypor;i++) {free(tabinter[i]);tabinter[i]=NULL;}
  free (tabinter);tabinter=NULL;
 }
 //printf("F5\n");
 //printf("hypoth= %d\n",hypoth);
 if (hypoth>0)
 {for (i=0;i<hypoth;i++) {free(tabhypo[i]);tabhypo[i]=NULL;}
  free (tabhypo);tabhypo=NULL;
 }
 //printf("F6\n");
 //printf("hypint= %d\n",hypint);
 if (hypint>0)
 {for (i=0;i<hypint;i++) {free(tabhypint[i]);tabhypint[i]=NULL;}
  free (tabhypint);tabhypint=NULL;
 }

//VG 09112006 Réinitialisation des variables à 0
/*
nbadd=0;
intercov=0;
nbhypor=0; 
nbhypor=0;
hypoth=0;
hypint=0;
printf("FIn\n");
*/

return(1);
}

void categorie()
{int i;
 nbcatego=1;
 suiv=base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {if (suiv->phen[0]>nbcatego) nbcatego=(int) suiv->phen[0];
   suiv=suiv->next;
 }

 nkat=nbcatego-1;
 nbsujktgo=(int *) malloc ((size_t) (nbcatego*sizeof(int)));
 for (i=0;i<nbcatego;i++)
 {nbsujktgo[i]=0;}

 suiv=base;
 if (msdata==0)
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {nbsujktgo[(int) suiv->phen[0]-1]+=(suiv->nblm==0);
  suiv=suiv->next;
 }
 else if (msdata==1)
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {nbsujktgo[(int) suiv->phen[0]-1]+=(suiv->nblm<nbloci-locmq+1);
  suiv=suiv->next;
 }
}

void polytomous()
{double val1,val2,val3,varest,det,res,denom,num,numm,probb;
 int hh1,hh2,v1,v2,nused=0,i,j,idx,k;
  char rep;
 double *dmat,probv0;


 dmat=(double *) malloc((size_t) (n*sizeof(double)));

 for (i=0;i<n;i++) mdvs2[i]=0;
 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}

 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {probv0=0;

  if (suiv->tnbhapo>0)
  {for (i=0;i<n;i++) {dmat[i]=0;}

   idx=((int) suiv->phen[0])-1;
   if (haplozero==0)
   {hh1=coding(suiv->hapest[0]);
    hh2=coding(suiv->hapest[1]);
   }

   denom=1;
   for (i=0;i<nkat;i++)
   {val2=2*effest[i];
    for (j=0;j<ajust;j++) val2+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];
   	if (haplozero==0)
	{if (hh1>0) {val2+=effest[hh1*nkat+i];}
     if (hh2>0) {val2+=effest[hh2*nkat+i];}
     for (j=0;j<nbadd;j++)
     {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
      {val2+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
     }
	 for (j=0;j<intercov;j++)
     {val2+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
     }

    }
    /*penser à nbadd+interecov*/
    denom+=exp(val2);
   }
   probv0=1/denom;

   if (idx==0)
   {for (i=0;i<nkat;i++)
    {res=2*effest[i];
     if (haplozero==0)
	 {if (hh1>0) res+=effest[hh1*nkat+i];
      if (hh2>0) res+=effest[hh2*nkat+i];
	  for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {res+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
      }
	  for (j=0;j<intercov;j++)
      {res+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) res+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];
     num=exp(res);
	 val1=-num*probv0;
     dmat[i]+=2*val1;
	 if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1*nkat+i]==1)) dmat[nitptp[hh1*nkat+i]]+=val1;
      if ((hh2>0) && (itptp[hh2*nkat+i]==1)) dmat[nitptp[hh2*nkat+i]]+=val1;

	  if (hypoth>0)
      {for (j=0;j<hypoth;j++)
       {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
        v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
        if ( ((hh1*nkat+i)==v2) && (nitptp[hh1*nkat+i]>-1)) {dmat[nitptp[hh1*nkat+i]]+=val1;}
	    if ( ((hh2*nkat+i)==v2) && (nitptp[hh2*nkat+i]>-1)) {dmat[nitptp[hh2*nkat+i]]+=val1;}
	   }
      }

	  if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][k][0],nkat).quot>0))
	     {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	       							           }
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	 								           }
	     }
	     else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
         {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                      }
          if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                      }
	     }
	     else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
         {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	     }
         else
	     {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][k][0]]]+=val1;
                                                dmat[nitptp[tabinter[j][0][0]]]-=val1;
		        						       }
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][k][0]]]+=val1;
  	                                            dmat[nitptp[tabinter[j][0][0]]]-=val1;
		  	 						           }
         }
        }
       }
      }





	  for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {dmat[nitptp[(nbhest+ajust)*nkat+j*nkat+i]]+=val1;}
      }
	  for (j=0;j<intercov;j++)
      {dmat[nitptp[(nbhest+ajust+nbadd)*nkat+j]]+=val1*suiv->z[tabint[j][1]-1]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) {dmat[nitptp[nbhest*nkat+j*nkat+i]]+=suiv->z[j]*val1;}
	}
   }
   else if (idx>0)
   {for (i=0;i<nkat;i++)
    {res=2*effest[i];
     if (haplozero==0)
	 {if (hh1>0) res+=effest[hh1*nkat+i];
      if (hh2>0) res+=effest[hh2*nkat+i];
      for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {res+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
      }
	  for (j=0;j<intercov;j++)
      {res+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) res+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];
     num=exp(res);
	 val1=-num*probv0;
     dmat[i]+=2*val1;
	 if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1*nkat+i]==1)) dmat[nitptp[hh1*nkat+i]]+=val1;
      if ((hh2>0) && (itptp[hh2*nkat+i]==1)) dmat[nitptp[hh2*nkat+i]]+=val1;

	  if (hypoth>0)
      {for (j=0;j<hypoth;j++)
       {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
        v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
        if ( ((hh1*nkat+i)==v2) && (nitptp[hh1*nkat+i]>-1)) {dmat[nitptp[hh1*nkat+i]]+=val1;}
	    if ( ((hh2*nkat+i)==v2) && (nitptp[hh2*nkat+i]>-1)) {dmat[nitptp[hh2*nkat+i]]+=val1;}
	   }
      }

	  if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][k][0],nkat).quot>0))
	     {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	       							           }
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	 								           }
	     }
	     else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
         {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                       }
          if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                       }
	     }
	     else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
         {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	     }
         else
	     {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][k][0]]]+=val1;
                                                dmat[nitptp[tabinter[j][0][0]]]-=val1;
		        						       }
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                                dmat[nitptp[tabinter[j][k][0]]]+=val1;
  	                                            dmat[nitptp[tabinter[j][0][0]]]-=val1;
		  	 						           }
         }
        }
       }
      }






	  for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {dmat[nitptp[(nbhest+ajust)*nkat+j*nkat+i]]+=val1;}
      }
	  for (j=0;j<intercov;j++)
      {dmat[nitptp[(nbhest+ajust+nbadd)*nkat+j]]+=val1*suiv->z[tabint[j][1]-1]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) {dmat[nitptp[nbhest*nkat+j*nkat+i]]+=suiv->z[j]*val1;}
	}
    dmat[idx-1]+=2;
    if (haplozero==0)
	{if ((hh1>0) && (itptp[hh1*nkat+idx-1]==1)) dmat[nitptp[hh1*nkat+idx-1]]+=1;
     if ((hh2>0) && (itptp[hh2*nkat+idx-1]==1)) dmat[nitptp[hh2*nkat+idx-1]]+=1;

     if (hypoth>0)
      {for (j=0;j<hypoth;j++)
       {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
        v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
        if ( ((hh1*nkat+idx-1)==v2) && (nitptp[hh1*nkat+idx-1]>-1)) {dmat[nitptp[hh1*nkat+idx-1]]+=1;}
	    if ( ((hh2*nkat+idx-1)==v2) && (nitptp[hh2*nkat+idx-1]>-1)) {dmat[nitptp[hh2*nkat+idx-1]]+=1;}
	   }
      }


	  if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][k][0],nkat).quot>0))
	     {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
	                                                dmat[nitptp[tabinter[j][k][0]]]+=1;
	       							               }
	      if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
	                                                dmat[nitptp[tabinter[j][k][0]]]+=1;
	 								               }
	     }
	     else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
         {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                   dmat[nitptp[tabinter[j][0][0]]]-=1;
	 	                                           }
          if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                    dmat[nitptp[tabinter[j][0][0]]]-=1;
	 	                                           }
	     }
	     else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
         {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;}
	      if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;}
	     }
         else
	     {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                    dmat[nitptp[tabinter[j][k][0]]]+=1;
                                                    dmat[nitptp[tabinter[j][0][0]]]-=1;
		      						          }
	      if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                    dmat[nitptp[tabinter[j][k][0]]]+=1;
  	                                                dmat[nitptp[tabinter[j][0][0]]]-=1;
		  			    				           }
         }
        }
       }
      }

	 for (j=0;j<nbadd;j++)
     {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
      {dmat[nitptp[(nbhest+ajust)*nkat+j*nkat+idx-1]]+=1;}
     }
	 for (j=0;j<intercov;j++)
     {dmat[nitptp[(nbhest+ajust+nbadd)*nkat+j]]+=suiv->z[tabint[j][1]-1]*(((hh1*nkat+idx-1)==tabint[j][0])+ ((hh2*nkat+idx-1)==tabint[j][0]));
     }
    }
    for (j=0;j<ajust;j++) {dmat[nitptp[nbhest*nkat+j*nkat+idx-1]]+=suiv->z[j];}
   }

   for (i=0;i<n;i++) {mdvs2[i]+=dmat[i];}
   for (i=0;i<n;i++) {for (j=0;j<n;j++) {mdvd2[i][j]+=dmat[i]*dmat[j];}}
  }
   suiv=suiv->next;
 }


 sysl(mdvd2,n);

 for (i=0;i<n;i++)  {modif2[i]=0;  for (j=0;j<n;j++) modif2[i]+=mdvd2[i][j]*mdvs2[j];}



 j=0;for (i=0;i<nall;i++) if (nitptp[i]>-1) {effest[i]+=modif2[nitptp[i]];j++;}


 if (hypoth>0)
 {for (j=0;j<hypoth;j++)
  {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
   v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
   if (div(v1,nkat).quot!=0) {effest[v2]=effest[v1];}
   else {effest[v2]=0;}
  }
 }


 if (interor==1)
  {for (j=0;j<nbhypor;j++)
   {for (i=1;i<nbor[j];i++)
    {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][i][0],nkat).quot>0))
	 {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]]+effest[tabinter[j][i][0]];
	 }
	 else if ( (div(tabinter[j][i][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
     {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]]-effest[tabinter[j][0][0]];
	 }
	 else if ( (div(tabinter[j][i][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
     {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]];
	 }
	 else
	 {effest[tabinter[j][i][1]]=effest[tabinter[j][0][1]]+effest[tabinter[j][i][0]]-effest[tabinter[j][0][0]];
	 }
 	}
   }
  }

 /* if (hypint>0)
 {for (i=0;i<hypint;i++)
  {effest[nbhest+ajust+nbadd+tabhypint[i][1]-1]=effest[nbhest+ajust+nbadd+tabhypint[i][0]-1];
  }
 }
*/
 for (i=0;i<n;i++) tempx+=(modif2[i]*modif2[i]);


 dmat=NULL; free (dmat);

}



/*Calcul de la matrice de variance/covaraince pour trait polytomique*/
void vpolyto(double **vintra)
{double val1,val2,val3,denom,res,num;
 int hh1,hh2,h1,h2,v1,v2,nused=0,i,j,idx,k;

 double *dmat,probv0;
 dindividu *parcour;


 dmat=(double *) malloc((size_t) ((n)*sizeof(double)));


 for (i=0;i<nall;i++) {for (j=0;j<nall;j++) mdvd2[i][j]=0;}

 suiv = base;
 while ((suiv!=NULL) && (suiv->next!=NULL))
 {probv0=0;

  if (suiv->tnbhapo>0)
  {for (i=0;i<n;i++) {dmat[i]=0;}

   idx=((int) suiv->phen[0])-1;
   if (haplozero==0)
   {hh1=coding(suiv->hapest[0]);
    hh2=coding(suiv->hapest[1]);
   }

   denom=1;
   for (i=0;i<nkat;i++)
   {val2=2*effest[i];
    for (j=0;j<ajust;j++) val2+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];

   	if (haplozero==0)
	{if (hh1>0) {val2+=effest[hh1*nkat+i];}
     if (hh2>0) {val2+=effest[hh2*nkat+i];}
     for (j=0;j<nbadd;j++)
     {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
      {val2+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
     }
	 for (j=0;j<intercov;j++)
     {val2+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
     }

    }
    denom+=exp(val2);
   }
   probv0=1/denom;

   if (idx==0)
   {for (i=0;i<nkat;i++)
    {res=2*effest[i];
     if (haplozero==0)
	 {if (hh1>0) res+=effest[hh1*nkat+i];
      if (hh2>0) res+=effest[hh2*nkat+i];
	  for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {res+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
      }
	  for (j=0;j<intercov;j++)
      {res+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) res+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];
     num=exp(res);
	 val1=-num*probv0;
     dmat[i]+=2*val1;
	 if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1*nkat+i]==1)) dmat[nitptp[hh1*nkat+i]]+=val1;
      if ((hh2>0) && (itptp[hh2*nkat+i]==1)) dmat[nitptp[hh2*nkat+i]]+=val1;

	  if (hypoth>0)
      {for (j=0;j<hypoth;j++)
       {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
        v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
        if ( ((hh1*nkat+i)==v2) && (nitptp[hh1*nkat+i]>-1)) {dmat[nitptp[hh1*nkat+i]]+=val1;}
	    if ( ((hh2*nkat+i)==v2) && (nitptp[hh2*nkat+i]>-1)) {dmat[nitptp[hh2*nkat+i]]+=val1;}
	   }
      }


	  if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][k][0],nkat).quot>0))
	     {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	       							           }
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	 								           }
	    }
	    else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
        {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                      }
         if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                      }
	    }
	    else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
        {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	     if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	    }
        else
	    {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][k][0]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
		      						          }
	     if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][k][0]]]+=val1;
  	                                           dmat[nitptp[tabinter[j][0][0]]]-=val1;
		  							          }
        }
       }
      }

     }





	  for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {dmat[nitptp[(nbhest+ajust)*nkat+j*nkat+i]]+=val1;}
      }
	  for (j=0;j<intercov;j++)
      {dmat[nitptp[(nbhest+ajust+nbadd)*nkat+j]]+=val1*suiv->z[tabint[j][1]-1]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) {dmat[nitptp[nbhest*nkat+j*nkat+i]]+=suiv->z[j]*val1;}
	}
   }
   else if (idx>0)
   {for (i=0;i<nkat;i++)
    {res=2*effest[i];
     if (haplozero==0)
	 {if (hh1>0) res+=effest[hh1*nkat+i];
      if (hh2>0) res+=effest[hh2*nkat+i];
      for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {res+=effest[(nbhest+ajust)*nkat+j*nkat+i];}
      }
	  for (j=0;j<intercov;j++)
      {res+=suiv->z[tabint[j][1]-1]*effest[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) res+=effest[nbhest*nkat+j*nkat+i]*suiv->z[j];
     num=exp(res);
	 val1=-num*probv0;
     dmat[i]+=2*val1;
	 if (haplozero==0)
     {if ((hh1>0) && (itptp[hh1*nkat+i]==1)) dmat[nitptp[hh1*nkat+i]]+=val1;
      if ((hh2>0) && (itptp[hh2*nkat+i]==1)) dmat[nitptp[hh2*nkat+i]]+=val1;

	  if (hypoth>0)
      {for (j=0;j<hypoth;j++)
       {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
        v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
        if ( ((hh1*nkat+i)==v2) && (nitptp[hh1*nkat+i]>-1)) {dmat[nitptp[hh1*nkat+i]]+=val1;}
	    if ( ((hh2*nkat+i)==v2) && (nitptp[hh2*nkat+i]>-1)) {dmat[nitptp[hh2*nkat+i]]+=val1;}
	   }
      }

     if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][k][0],nkat).quot>0))
	     {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	       							           }
	      if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
	                                            dmat[nitptp[tabinter[j][k][0]]]+=val1;
	 								           }
	    }
	    else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
        {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                      }
         if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
	 	                                      }
	    }
	    else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
        {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	     if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;}
	    }
        else
	    {if ((hh1*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][k][0]]]+=val1;
                                               dmat[nitptp[tabinter[j][0][0]]]-=val1;
		      						          }
	     if ((hh2*nkat+i)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=val1;
                                               dmat[nitptp[tabinter[j][k][0]]]+=val1;
  	                                           dmat[nitptp[tabinter[j][0][0]]]-=val1;
		  							          }
         }
        }
       }
      }
	  for (j=0;j<nbadd;j++)
      {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
       {dmat[nitptp[(nbhest+ajust)*nkat+j*nkat+i]]+=val1;}
      }
	  for (j=0;j<intercov;j++)
      {dmat[nitptp[(nbhest+ajust+nbadd)*nkat+j]]+=val1*suiv->z[tabint[j][1]-1]*(((hh1*nkat+i)==tabint[j][0])+ ((hh2*nkat+i)==tabint[j][0]));
      }
	 }
	 for (j=0;j<ajust;j++) {dmat[nitptp[nbhest*nkat+j*nkat+i]]+=suiv->z[j]*val1;}
	}
    dmat[idx-1]+=2;
    if (haplozero==0)
	{if ((hh1>0) && (itptp[hh1*nkat+idx-1]==1)) dmat[nitptp[hh1*nkat+idx-1]]+=1;
     if ((hh2>0) && (itptp[hh2*nkat+idx-1]==1)) dmat[nitptp[hh2*nkat+idx-1]]+=1;

      if (hypoth>0)
      {for (j=0;j<hypoth;j++)
       {v1=tabhypo[j][0]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][1]*(tabhypo[j][1]<tabhypo[j][0]);
        v2=tabhypo[j][1]*(tabhypo[j][0]<tabhypo[j][1])+tabhypo[j][0]*(tabhypo[j][1]<tabhypo[j][0]);
        if ( ((hh1*nkat+idx-1)==v2) && (nitptp[hh1*nkat+idx-1]>-1)) {dmat[nitptp[hh1*nkat+idx-1]]+=1;}
	    if ( ((hh2*nkat+idx-1)==v2) && (nitptp[hh2*nkat+idx-1]>-1)) {dmat[nitptp[hh2*nkat+idx-1]]+=1;}
	   }
      }

      if (interor==1)
      {for (j=0;j<nbhypor;j++)
       {for (k=1;k<nbor[j];k++)
        {if ( (div(tabinter[j][0][0],nkat).quot==0) && (div(tabinter[j][k][0],nkat).quot>0))
	     {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
	                                                dmat[nitptp[tabinter[j][k][0]]]+=1;
	       							               }
	      if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
	                                                dmat[nitptp[tabinter[j][k][0]]]+=1;
	 								               }
	    }
	    else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot>0))
        {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                   dmat[nitptp[tabinter[j][0][0]]]-=1;
	 	                                          }
         if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                   dmat[nitptp[tabinter[j][0][0]]]-=1;
	 	                                          }
	    }
	    else if ( (div(tabinter[j][k][0],nkat).quot==0) && (div(tabinter[j][0][0],nkat).quot==0))
        {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;}
	     if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;}
	    }
        else
	    {if ((hh1*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                   dmat[nitptp[tabinter[j][k][0]]]+=1;
                                                   dmat[nitptp[tabinter[j][0][0]]]-=1;
		      						              }
	     if ((hh2*nkat+idx-1)==tabinter[j][k][1]) {dmat[nitptp[tabinter[j][0][1]]]+=1;
                                                   dmat[nitptp[tabinter[j][k][0]]]+=1;
  	                                               dmat[nitptp[tabinter[j][0][0]]]-=1;
		  							              }
         }
        }
       }
      }


	 for (j=0;j<nbadd;j++)
     {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
      {dmat[nitptp[(nbhest+ajust)*nkat+j*nkat+idx-1]]+=1;}
     }
	 for (j=0;j<intercov;j++)
     {dmat[nitptp[(nbhest+ajust+nbadd)*nkat+j]]+=suiv->z[tabint[j][1]-1]*(((hh1*nkat+idx-1)==tabint[j][0])+ ((hh2*nkat+idx-1)==tabint[j][0]));
     }
    }
    for (j=0;j<ajust;j++) {dmat[nitptp[nbhest*nkat+j*nkat+idx-1]]+=suiv->z[j];}
   }
   for (i=0;i<n;i++) {for (j=0;j<n;j++) {mdvd2[i][j]+=dmat[i]*dmat[j];}}
  }
   suiv=suiv->next;
 }

 sysl(mdvd2,n);
 for (i=0;i<nall;i++)
	for (j=0;j<nall;j++)
	  if ((nitptp[i]>-1) && (nitptp[j]>-1))  vintra[i][j]+=mdvd2[nitptp[i]][nitptp[j]];

 dmat=NULL; free (dmat);

}



/*********************CALCUL DE LA VRAISEMBLANCE TOTALE ******************************************/
double likepoly(double *frqsem, double *lgoddsem)
{double val1,v1,vrais,h1,h2,like,p1,vraisfa=0,vv2,denom,num,val2;
 int hh1,hh2/*,idx*/,i,j,ii,idx;
  suiv = base;

  while ((suiv!=NULL) && (suiv->next!=NULL))       /*ATTENTION AU PB EVENTUEL DES DONNES MANQUANTES*/
  {val1=suiv->phen[0];
   if (suiv->tnbhapo>0)
   {vrais=0;
    idx=((int) suiv->phen[0])-1;
	for (i=0;i<suiv->tnbhapo;i++)
    {h1=frqsem[suiv->idnb[i][0]];h2=frqsem[suiv->idnb[i][1]];
     if ( (h1>0) && (h2>0))
     {p1=h1*h2*(2-(suiv->idnb[i][0]==suiv->idnb[i][1]));
   	  if (haplozero==0)
      {hh1=coding(suiv->idnb[i][0]);
       hh2=coding(suiv->idnb[i][1]);
      }
      denom=1;
      for (ii=0;ii<nkat;ii++)
	  {val2=2*lgoddsem[ii];
       for (j=0;j<ajust;j++) val2+=lgoddsem[nbhest*nkat+j*nkat+ii]*suiv->z[j];
       if (haplozero==0)
	   {if (hh1>0) {val2+=lgoddsem[hh1*nkat+ii];}
        if (hh2>0) {val2+=lgoddsem[hh2*nkat+ii];}
        for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {val2+=lgoddsem[(nbhest+ajust)*nkat+j*nkat+ii];}
        }
		for (j=0;j<intercov;j++)
        {val2+=suiv->z[tabint[j][1]-1]*lgoddsem[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+ii)==tabint[j][0])+ ((hh2*nkat+ii)==tabint[j][0]));
        }
	   }
       denom+=exp(val2);
      }
      num=1;
      if (idx>0)
      {val2=2*lgoddsem[idx-1];
       for (j=0;j<ajust;j++) val2+=lgoddsem[nbhest*nkat+j*nkat+idx-1]*suiv->z[j];
       if (haplozero==0)
       {if (hh1>0) {val2+=lgoddsem[nkat*hh1+idx-1];}
        if (hh2>0) {val2+=lgoddsem[nkat*hh2+idx-1];}
        for (j=0;j<nbadd;j++)
        {if ( ((tadd[j][0]-1==hh1) && (tadd[j][1]-1==hh2)) || ((tadd[j][0]-1==hh2) && (tadd[j][1]-1==hh1)) )
         {val2+=lgoddsem[(nbhest+ajust)*nkat+j*nkat+idx-1];}
        }
	    for (j=0;j<intercov;j++)
        {val2+=suiv->z[tabint[j][1]-1]*lgoddsem[(nbhest+ajust+nbadd)*nkat+j]*(((hh1*nkat+idx-1)==tabint[j][0])+ ((hh2*nkat+idx-1)==tabint[j][0]));
        }
	   }
       /*Penser à nbadd+intercov*/
        num=exp(val2);
      }
      like=num/denom;
      vrais+=p1*like;
	 }
    }
    if (vrais>0) {vraisfa -= log(vrais);}
   }
   suiv= suiv->next;
  }
  return(-vraisfa);
}




