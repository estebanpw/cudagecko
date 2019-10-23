/* lista.c

functions to handle a list



*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "lista.h"

void Insertar(Lista *lista, int v,Fragmentv3 f,int n_event) {
   pNodo nuevo, actual;

   /* Crear un nodo nuevo */
   nuevo = (pNodo)malloc(sizeof(tipoNodo));
   nuevo->valor = v;
   nuevo->f=f;
   nuevo->n_event=n_event;
   
   /* Colocamos actual en la primera posición de la lista */
   actual = *lista;
   if(actual) while(actual->anterior) actual = actual->anterior;
   
   /* Si la lista está vacía o el primer miembro es mayor que el nuevo */
   if(!actual || actual->valor > v) {
      /* Añadimos la lista a continuación del nuevo nodo */
      nuevo->siguiente = actual; 
      nuevo->anterior = NULL;
      if(actual) actual->anterior = nuevo;
      if(!*lista) *lista = nuevo;
   }
   else {
      /* Avanzamos hasta el último elemento o hasta que el siguiente tenga 
         un valor mayor que v */
      while(actual->siguiente && actual->siguiente->valor <= v) 
         actual = actual->siguiente;
      /* Insertamos el nuevo nodo después del nodo anterior */
      nuevo->siguiente = actual->siguiente;
      actual->siguiente = nuevo;
      nuevo->anterior = actual;
      if(nuevo->siguiente) nuevo->siguiente->anterior = nuevo;
   }
}

void Borrar(Lista *lista, int v) {
   pNodo nodo;
   
   
   /* Buscar el nodo de valor v */
   nodo = *lista;
   while(nodo && nodo->valor <v) nodo = nodo->siguiente;
   while(nodo && nodo->valor > v) nodo = nodo->anterior;

   /* El valor v no está en la lista */
   if(!nodo || nodo->valor != v) return;
   
   /* Borrar el nodo */
   /* Si lista apunta al nodo que queremos borrar, apuntar a otro */
   if(nodo == *lista){
     if(nodo->anterior) {*lista = nodo->anterior;}
     else {*lista = nodo->siguiente;}
   }
   if(nodo->anterior) /* no es el primer elemento */
      nodo->anterior->siguiente = nodo->siguiente;
   if(nodo->siguiente) /* no es el último nodo */
      nodo->siguiente->anterior = nodo->anterior;
   free(nodo);
}
void BorrarLista(Lista *lista) {
   pNodo nodo, actual;

   actual = *lista;
   while(actual->anterior) actual = actual->anterior;

   while(actual) {
      nodo = actual;
      actual = actual->siguiente;
      free(nodo);
   }
   *lista = NULL;
}

void MostrarLista(Lista lista, int orden,char* nombre,int num) {
   pNodo nodo = lista;
	int i=0;
   if(!lista) printf("Lista vacía");

   nodo = lista;
   if(orden == 1) {
      while(nodo->anterior) nodo = nodo->anterior;
     // printf("Orden ascendente: ");
      while(nodo) {
        // printf("%d -> ", nodo->valor);
		printFrag(nodo->f,nombre,i++,num,nodo->n_event);
         nodo = nodo->siguiente;
      }
   }
   else {
      while(nodo->siguiente) nodo = nodo->siguiente;
     // printf("Orden descendente: ");
      while(nodo) {
        printFrag(nodo->f,nombre,i++,num,nodo->n_event);
         nodo = nodo->anterior;
      }
   }
   
   //printf("\n");
}

void printFrag(Fragmentv3 f,char* nombre,int i,int num,int n_event){

	
	fprintf(NULL,"%d\t%s\t%d\t%ld\t%ld\t%ld\t%ld\t%f\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\n",i,nombre,(int)num,f.xIni,f.xFin,f.yIni,f.yFin,f.similarity,(int)f.id,(int)f.seqY,(int)f.seqX,(int)f.block,(int)f.gap,f.strand,(int)f.reversions,(int)f.translocations,(int)f.duplications,(int)f.events,(int)f.score);
	printf("%ld\t%ld\t%d\n",f.seqY,f.seqX,n_event);
		

}


void CopiarListas(Lista* origen,Lista* dest){ // SUponemos las dos listas mismo tamaño

	pNodo onodo = *origen;
	pNodo dnodo = *dest;
	while(onodo->anterior && dnodo->anterior) {onodo = onodo->anterior;dnodo = dnodo->anterior;}

   while(onodo && dnodo) {
   /*
	dnodo->f=*(onodo->f);
	dnodo->valor=*(onodo->valor);
      */
	  
	  dnodo = dnodo->siguiente;
      onodo = onodo->siguiente;
     
   }
	

}

