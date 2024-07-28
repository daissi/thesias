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

#include "thesiaslib.h"

extern double igam ( double, double );
extern double igami ( double, double );
extern double polevl ( double, double *, int );
extern double ndtri (double);
extern double gamma (double);
extern double chdtrc(double, double);
extern int mtherr(char *, int );
extern int thesiasRun(char*, int, int, int*, int, int, int, int, int, int, int, int, int, int*, int, int);

JNIEXPORT jint JNICALL Java_thesiaslib_thesiasRun
  (JNIEnv      *env, 
   jobject     obj, 
   jstring     sFileName, 
   jint        maxvarfic,
   jint        nbloci,
   jintArray   idloci,
   jint        ldmatrix,
   jint        msdata,
   jint        R2,
   jint        chxt,
   jint        num0,
   jint        idtime,
   jint        offset,
   jint        idoffset,
   jint        ajust,
   jintArray   numajust,
   jint	       xlnk,
   jint        numsx)

{
	int i;
	jint *aNdloci = (*env)->GetIntArrayElements(env, idloci, 0);
	jint *aNumajust = (*env)->GetIntArrayElements(env, numajust, 0);
	const char *cFileName = (*env)->GetStringUTFChars(env, sFileName, 0);

	printf("Entrée dans le pivot JNI: appel natif de thesiasRun VGT \n");
	i = thesiasRun(cFileName,
		       maxvarfic,
		       nbloci,
		       aNdloci,
		       ldmatrix,
		       msdata,
		       R2,
		       chxt,
		       num0,
		       idtime,
		       offset,
		       idoffset,
		       ajust,
		       aNumajust,
		       xlnk,
		       numsx);

	(*env)->ReleaseIntArrayElements(env, idloci, aNdloci, 0);
	(*env)->ReleaseIntArrayElements(env, numajust, aNumajust, 0);
	(*env)->ReleaseStringUTFChars(env, sFileName, cFileName);

	return i;
}


JNIEXPORT jint JNICALL
JNI_OnLoad(JavaVM *vm, void *reserved){

	return JNI_VERSION_1_4;

}

void JNICALL
JNI_OnUnload(JavaVM *vm, void *reserved){
	return;
}


