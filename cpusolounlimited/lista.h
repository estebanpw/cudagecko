/*lista.h


*/


#ifndef LISTA_H_
#define LISTA_H_

#include "fragmentv3.h"

typedef struct _nodo {
	Fragmentv3 f;
   int valor;
   int n_event;
   struct _nodo *siguiente;
   struct _nodo *anterior;
} tipoNodo;
 
typedef tipoNodo *pNodo;
typedef tipoNodo *Lista;


void Insertar(Lista *lista, int v,Fragmentv3 f,int n_event);
void Borrar(Lista *lista, int v);
void printFrag(Fragmentv3 f,char* nombre,int i,int num,int n_event);
void BorrarLista(Lista *lista);
void CopiarListas(Lista* origen,Lista* dest); // SUponemos las dos listas mismo tamaño




#endif /* LISTA_H_ */
