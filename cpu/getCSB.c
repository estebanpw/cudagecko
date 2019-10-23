/* getCSB.c



*/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "fragmentv3.h"
#include "fragmentv2.h"
#include "lista.h"
#define  MAX_LEVELS  900000

/* Subprograms */

void quickSortX(Fragmentv3 *f, int elements);
void quickSortScore(Fragmentv3 *f, int elements);
void quickSortY(Fragmentv3 *f, int elements, int ytotal);
void reorder( Fragmentv3* f, int nf, int ytotal);

void resetFrags(Fragmentv3* f,int nf);

Lista ordenarXY(Lista *lista,char c);
void ordenarLista(Lista *lista);

int S(Fragmentv3 f); // It returns '1' for 'f' strand and '-1' for 'r' strand;
int C(Fragmentv3 a, Fragmentv3 ci,Fragmentv3 cd, Fragmentv3 s); // Collinearity of ant, current, sig.
int L(Fragmentv3 a, Fragmentv3 c); // Linear of ant and current
int CSB(Lista *lista,Fragmentv3 *frags); // Join all posible CSB. return the number of CSB

int NewList(Lista* lista,Fragmentv3* f,int nf);
void saveCSB(Lista lista,int xtotal,int ytotal);
void saveFragments(Fragmentv3* f,int nf,int xtotal,int ytotal);

/********************/

/* Globals */
//int dist;
//int size;
char** AV;

//int MOSTRAR=1;
//int CONTADOR=0;

// Constantes extensión de gap
int INI=0;
float EXT=3;




int main(int ac,char** av){

	AV=av;
	if(ac<4){
		printf("Use: getCSB fragsv3.fil csb.fv3 csb.txt \n");
		exit(-1);
	}
	/***********************/
	/* General vars */
	int size_of_list;
	//int i;


	/* Read Fragments */
		Fragmentv3* f;
		int nf,xtotal,ytotal;
		f=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);

		
	/* Reset frags */	
	resetFrags(f,nf);
	
	/* Order frags */
	reorder(f,nf,ytotal); 
	
	/* Creates a list */
	Lista lista = NULL;
	size_of_list=NewList(&lista,f,nf);
//+++++++++++++++++++++++++++	
/************ START PROGRAM ***************************************************************/
	int n=1;
	
		
	// Get CSBs
	if(size_of_list>1){
		n=CSB(&lista,f);
	}
	
	
    saveFragments(f,nf,xtotal,ytotal);

return 0;	
}








/*************************/
int NewList(Lista* lista,Fragmentv3* f,int nf){
		int i;
		int size_of_list=0;
		int m;
		m=0;
		for(i=0;i<nf;i++){
			
			if(f[i].block>=0){
				size_of_list++;
				
				f[i].block=i+1;
				f[i].id=i;
				Insertar(lista,m++,f[i],0);
			
			}else{
				f[i].id=-1;
			}
			//printf("%d\t%ld\t%ld\t%ld\t%d\n",i,f[i].id,f[i].score,f[i].length,f[i].block);
			
		}
		return size_of_list;
	}
/***************************************************/
/***************************************************/
void saveFragments(Fragmentv3 *frags,int nf,int xtotal,int ytotal){

FILE* fs;
	if((fs=fopen(AV[2],"wb"))==NULL){
		printf("***ERROR Opening output file %s\n",AV[2]);
		exit(-1);
	}
FILE* ftxt;
	if((ftxt=fopen(AV[3],"w"))==NULL){
		printf("***ERROR Opening output file %s\n",AV[3]);
		exit(-1);
	}
//printf("abierto %s\n",av[4]);
Fragmentv3 fragv3;

	uint64_t xtotal64,ytotal64;
	xtotal64 = xtotal;
	ytotal64 = ytotal;
	fwrite(&xtotal,sizeof(int),1,fs);
	fwrite(&ytotal,sizeof(int),1,fs);
//	fwrite(&xtotal,sizeof(int),1,fs);
//	fwrite(&ytotal,sizeof(int),1,fs);



	
	int i;
	for(i=0;i<nf;i++){
		fragv3=frags[i];
		fwrite(&fragv3,sizeof(Fragmentv3),1,fs);
		fprintf(ftxt,"%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%d\n",fragv3.xIni,fragv3.yIni,fragv3.xFin,fragv3.yFin,fragv3.length,fragv3.strand,fragv3.score,fragv3.block);

		
	}
	
	fclose(fs);
	fclose(ftxt);
	
	//getchar();
}
/****************/
void saveCSB(Lista lista,int xtotal,int ytotal){

FILE* fs;
	if((fs=fopen(AV[2],"wb"))==NULL){
		printf("***ERROR Opening output file %s\n",AV[2]);
		exit(-1);
	}
FILE* ftxt;
	if((ftxt=fopen(AV[3],"w"))==NULL){
		printf("***ERROR Opening output file %s\n",AV[3]);
		exit(-1);
	}
//printf("abierto %s\n",av[4]);
Fragmentv3 fragv3;

	uint64_t xtotal64,ytotal64;
	xtotal64 = xtotal;
	ytotal64 = ytotal;
	fwrite(&xtotal,sizeof(int),1,fs);
	fwrite(&ytotal,sizeof(int),1,fs);
//	fwrite(&xtotal,sizeof(int),1,fs);
//	fwrite(&ytotal,sizeof(int),1,fs);


	pNodo nodo;
	nodo=lista;
	printf("-------------------------\n");
	fprintf(ftxt,"xIni\tyIni\txFin\tyFin\tlength\tstrand\tscore\n");
	
	while(nodo && nodo->siguiente){
		fragv3 = nodo->f;
		fwrite(&fragv3,sizeof(Fragmentv3),1,fs);
		fprintf(ftxt,"%ld\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%d\n",fragv3.xIni,fragv3.yIni,fragv3.xFin,fragv3.yFin,fragv3.length,fragv3.strand,fragv3.score,fragv3.block);

		nodo= nodo->siguiente;
	}
	
	fclose(fs);
	fclose(ftxt);
	printf("------------FINJ-------------\n");
	//getchar();
}
/***************************************************/
int S(Fragmentv3 f){

	if(f.strand=='f')return 1;
	if(f.strand=='r')return -1;
	return 0;
}
/***************************************************/
int C(Fragmentv3 a, Fragmentv3 ci, Fragmentv3 cd, Fragmentv3 s){

	if( L(a,ci) && L(cd,s) )return 1;
	return 0;
}
/***************************************************/
int L(Fragmentv3 a, Fragmentv3 c){

	int A,C;
	A=a.seqY;
	C=c.seqY;
	if( ( (A-S(a)) == (C-2*S(c)) ) && (S(c) == S(a))  )return 1;
	return 0;
}
/***************************************************/
int CSB(Lista *lista,Fragmentv3 *frags){
	int n=1;


	int penal=0;
	int joined=1;
	
	pNodo nodo; 
	while(joined){
		joined = 0;
		nodo = *lista;
		while(!nodo->anterior)nodo=nodo->siguiente;
		while(nodo->siguiente) {
		
		
		
				if( L(nodo->anterior->f,nodo->f) ){
					
					penal = INI + EXT*(abs( nodo->anterior->f.xFin- nodo->f.xIni)+abs( nodo->anterior->f.yFin- nodo->f.yIni)); 
					//printf("penal:%d\n",penal);
					
					nodo->anterior->f.xFin = nodo->f.xFin;
					nodo->anterior->f.yFin = nodo->f.yFin;
					//nodo->anterior->f.score = nodo->anterior->f.score + nodo->f.score;
					nodo->anterior->f.score = nodo->anterior->f.score + nodo->f.score- penal;
					nodo->anterior->f.length = abs(nodo->anterior->f.yFin - nodo->anterior->f.yIni);
					
					nodo->anterior->f.seqY = nodo->f.seqY;
					nodo->anterior->f.seqX = nodo->f.seqX;
					//nodo->anterior->f.block=n;
					//frags[nodo->f.id].block=n;
					//frags[nodo->anterior->f.id].block=n;
					frags[nodo->f.id].block=frags[nodo->anterior->f.id].block;
					
					
					
					joined = 1;
					Borrar(lista, nodo->valor);
					
				}else{
				
					
					nodo = nodo->siguiente;
					n++;
					
				}
		}
		// Para el ultimo
		if( L(nodo->anterior->f,nodo->f) ){
			penal = INI + EXT*(abs( nodo->anterior->f.xFin- nodo->f.xIni)+abs( nodo->anterior->f.yFin- nodo->f.yIni)); 
			//printf("penal:%d\n",penal);
			
			nodo->anterior->f.xFin = nodo->f.xFin;
			nodo->anterior->f.yFin = nodo->f.yFin;
			//nodo->anterior->f.score = nodo->anterior->f.score + nodo->f.score;
			nodo->anterior->f.score = nodo->anterior->f.score + nodo->f.score- penal;
			nodo->anterior->f.length = abs(nodo->anterior->f.yFin - nodo->anterior->f.yIni);
			
			nodo->anterior->f.seqY = nodo->f.seqY;
			nodo->anterior->f.seqX = nodo->f.seqX;
			
			frags[nodo->f.id].block=frags[nodo->anterior->f.id].block;
					
			joined = 1;	
			Borrar(lista, nodo->valor);
		}
		if(n>=1)ordenarLista(lista);
	}
	return n;
}
/***************************************************/
Lista ordenarXY(Lista *lista,char c){

	Lista nueva=NULL;
	pNodo nodo= *lista;
	nodo = *lista;
/*******************************************/	
/* Primero buscamos el valor menor, para que sea el primero en insertarse */

	int valor=0;
	pNodo primero=nodo;
	
	//while(nodo->anterior)nodo=nodo->anterior;
	if(c=='x'){valor=nodo->f.seqX;primero=nodo;}
	if(c=='y'){valor=nodo->f.seqY;primero=nodo;}
	
	
	while(nodo){
		if(c=='x'){
			if(valor>=nodo->f.seqX && nodo->f.block>=0){
				valor=nodo->f.seqX;
				primero=nodo;
			}
		
		}else{
			if(valor>=nodo->f.seqY && nodo->f.block>=0){
				valor=nodo->f.seqY;
				primero=nodo;
			}
		
		}
		nodo=nodo->siguiente;
	}
//	if(c=='x'){printf("aqui primero->f.seqX: %d valor:%d\n",primero->f.seqX,valor);}
//	if(c=='y'){printf("aqui primero->f.seqY: %d valor:%d\n",primero->f.seqY,valor);}
	

	//Insertamos el primero 
	if(c=='x'){
//		printf("X insertamos en la posicion %d\n",primero->f.seqX);
		Insertar(&nueva,primero->f.seqX,primero->f,primero->n_event);
		Borrar(lista,primero->valor);
	}else{
///		printf("Y insertamos en la posicion %d\n",primero->f.seqY);
		Insertar(&nueva,primero->f.seqY,primero->f,primero->n_event);
//		printf("ya hemos insertado\n");
		Borrar(lista,primero->valor);
	}

/*******************************************/
	
	nodo = *lista;
	while(nodo->anterior)nodo=nodo->anterior;
	while(nodo){
		// Insertarmos en el orden de XY
		if(c=='x'){
			Insertar(&nueva,nodo->f.seqX,nodo->f,nodo->n_event);
		}else{
			Insertar(&nueva,nodo->f.seqY,nodo->f,nodo->n_event);
			
		}
		nodo = nodo->siguiente;
	}

	
	pNodo newNodo = nueva;
	newNodo = nueva;
	int i=0;
	while(newNodo){
		if(c=='x'){
			
			newNodo->f.seqX=i++;
		}else{
			newNodo->f.seqY=i++;
		}
		
		newNodo = newNodo->siguiente;
	}

	return nueva;

	
}
/****************************************/
void reorder( Fragmentv3* f, int nf, int ytotal){
	int i;
	int j;

	/* Sort Y*/
	quickSortY(&f[0],nf,ytotal);
	j=0;
	for(i=0;i<nf;i++){
		if(f[i].block>=0)f[i].seqY=j++;
	}
		

	/* Sort X*/
	quickSortX(&f[0],nf);
	j=0;
	for(i=0;i<nf;i++){
		if(f[i].block>=0)f[i].seqX=j++;
	}
	

}
/******/
void resetFrags(Fragmentv3* f, int nf){
int i;
for(i=0;i<nf;i++){
	f[i].gap=0;

	f[i].events=0;
	f[i].id=i;
	
	f[i].reversions=0;

}

}
/**** Ordenar lista ****/
void ordenarLista(Lista *lista){

			// Ordenamos
			Lista nueva = NULL;
			nueva=ordenarXY(lista,'y');
			BorrarLista(lista);
			*lista=nueva;
			
			
			
			// Ahora por X
			Lista nuevaX = NULL;
			nuevaX=ordenarXY(lista,'x');
			BorrarLista(lista);
			*lista=nuevaX;
}

/****************************************/

//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7
void quickSortScore(Fragmentv3 *f, int elements) {

  


  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].score>=piv.score && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].score<=piv.score && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}}
/*****************/
void quickSortX(Fragmentv3 *f, int elements) {

  


  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].xIni>=piv.xIni && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].xIni<=piv.xIni && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}}
/*****************/
void quickSortY(Fragmentv3 *f, int elements,int ytotal) {

	
	
  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
	
/*	
// Change Y component
long ytemp;
	for(i=0;i<elements;i++){
		if(f[i].strand=='r'){
			//f[i].yIni=ytotal-f[i].yIni;
			//f[i].yFin=ytotal-f[i].yFin;

			ytemp=f[i].yIni;
			f[i].yIni=f[i].yFin;
			f[i].yFin=ytemp;	
		}
	}

*/
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].yIni>=piv.yIni && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].yIni<=piv.yIni && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}

// Change Y component
/*
	for(i=0;i<elements;i++){
		if(f[i].strand=='r'){
			ytemp=f[i].yIni;
			f[i].yIni=f[i].yFin;
			f[i].yFin=ytemp;	

			//f[i].yIni=ytotal-f[i].yIni;
			//f[i].yFin=ytotal-f[i].yFin;


		}
	}
*/

}

