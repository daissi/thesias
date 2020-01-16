				The THESIAS program

			"Testing Haplotype Effects In Association Studies"



		by Tregouet David-Alexandre (david.tregouet@chups.jussieu.fr)
				INSERM U525,
	      Genetic Epidemiology and Molecular of Cardiovascular Pathologies
				Paris, France


/*THESIAS program was updated on Tuesday 21st february 2006 :
Haplotype analysis of a categorical phenotype is now possible by use of a polytomous formulation.
Besides, X-linked SNPs analysis is now available for case-control and quantitative analysis.

*/


**************************************************************************************************
The objectif of the THESIAS program is to performed haplotype-based association analysis 
in unrelated individuals. This program is based on the maximum likelihood model described
in Tregouet et al. 2002 and is linked to the SEM algorithm (Tregouet et al. 2004).
THESIAS allows the simultaneous estimation of haplotype frequencies and of their associated
effects on the phenotype of interest. In its current version, quantitative, qualitative, categorical
and survival outcome can be studied.
Covariate-adjusted haplotype effects as well as haplotype x covariate interactions can 
be investigated. 


Every information on how to run THESIAS and on what THESIAS do can be found in the Documentation3.1.htm File.

In order to illustrate the use of the THESIAS program, several .txt and .html files are available
for download through the example.zip file. 

**************************************************************************************************
The Data.txt file includes one quantitative phenotype (in position 2),two polymorphisms 
(in position 3 and 4) and two additional covariates (in position 5 and 6).

In the Answer0.txt file are given the appropriate responses to be given to THESIAS in order to get 
the parameter file (para.txt). The results files should then be similar to para0.txt, result0.txt 
and result0.htm 

In the Answer1.txt file are given the appropriate responses to be given to THESIAS in order to get
the loglikelihood of a model testing the effect of the two covariates on the phenotype but with no 
haplotype effects. The results files should then be similar to result1.txt and result1.htm .


If one wants to further estimate haplotype effects, para.txt must have to be modified as in para1.txt .
Then run the THESIAS program with the responses given in the answer2.txt file.  The results files should
then be similar to result2.txt and result2.htm .

**************************************************************************************************

In order to illustrate the use of the extended THESIAS program for survival outcome analysis,
several .txt files are available for download through the example2.zip file. 

The data.dat file includes one censored phenotype (in position 3) with its associated survival time 
(in position 4),six polymorphisms (in position 5 to 10) and three additional covariates
(in position 11 to 13).

In the Answer0.txt file are given the appropriate responses to be given to THESIAS in order to get 
the parameter file (para.txt). The results files should then be similar to para0.txt, result0.txt 

In the Answer1.txt file are given the appropriate responses to be given to THESIAS in order to get
the loglikelihood of a model testing the effects of the five most frequent haplotypes. The results file
should then be similar to result1.txt.

**************************************************************************************************


 

 