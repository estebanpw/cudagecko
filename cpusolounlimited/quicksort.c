#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksort.h"

#define BUF_SIZE 1024
#define SWAP(a,b,t) t=a; a=b; b=t;
#define STRING_SIZE 1024

typedef struct {
	BaseType* a;
	int l;
	int r;
	int nth;
} PQSortArgs;

typedef struct {
	char fin1[STRING_SIZE];
	char fin2[STRING_SIZE];
	char fout[STRING_SIZE];
} PMergeArgs;

void assertNotNull(void* p, char* msg){
	if(p==NULL){
		terror(msg);
	}
}

void assertIntEQ(int a, int b, char* msg){
	if(a!=b){
		terror(msg);
	}
}

void assertIntGE(int a, int b, char* msg){
	if(a<b){
		terror(msg);
	}
}

int bufMerge(BaseType* a1, int n1, BaseType* a2, int n2, BaseType* m){
	int i1=0,i2=0;
	int j=0;

	while(i1<n1 && i2<n2){
		if(GT(a1[i1],a2[i2])){
			m[j]=a2[i2];
			i2++;
			j++;
		}else{
			m[j]=a1[i1];
			i1++;
			j++;
		}
	}

	while(i1<n1){
		m[j]=a1[i1];
		i1++;
		j++;
	}

	while(i2<n2){
		m[j]=a2[i2];
		i2++;
		j++;
	}

	return j;
}

void *PMerge(void* a){
	PMergeArgs* args = (PMergeArgs*)a;

	FILE* fin1 = fopen(args->fin1,"rb");
	assertNotNull(fin1,args->fin1);

	FILE* fin2 = fopen(args->fin2,"rb");
	assertNotNull(fin2,args->fin2);

	FILE* fout = fopen(args->fout,"wb");
	assertNotNull(fout,args->fout);

	BaseType *a1 = (BaseType*)calloc(BUF_SIZE,sizeof(BaseType));
	assertNotNull(a1,"calloc");

	BaseType *a2 = (BaseType*)calloc(BUF_SIZE,sizeof(BaseType));
	assertNotNull(a2,"calloc");

	BaseType *m = (BaseType*)calloc(2*BUF_SIZE,sizeof(BaseType));
	assertNotNull(m,"calloc");

	int i,j,k;
	int r1,r2,tmp;

	r1=fread(a1,sizeof(BaseType),BUF_SIZE,fin1);
	r2=fread(a2,sizeof(BaseType),BUF_SIZE,fin2);

	i=j=k=0;
	while(r1>0 && r2>0){
		while(i<r1 && j<r2){
			if(GT(a1[i],a2[j])){
				m[k]=a2[j];
				j++;
			}else{
				m[k]=a1[i];
				i++;
			}
			k++;
			if(k>=2*BUF_SIZE){
				tmp=fwrite(m,sizeof(BaseType),k,fout);
		    assertIntEQ(tmp,k,"fwrite");
				k=0;
			}
		}

		if(i>=r1){
			r1=fread(a1,sizeof(BaseType),BUF_SIZE,fin1);
			i=0;
		}

		if(j>=r2){
			r2=fread(a2,sizeof(BaseType),BUF_SIZE,fin2);
			j=0;
		}
	}

	if(k>0){
		tmp=fwrite(m,sizeof(BaseType),k,fout);
		assertIntEQ(tmp,k,"fwrite");
	}

	if(i<r1){
		tmp=fwrite(a1+i,sizeof(BaseType),r1-i,fout);
		assertIntEQ(tmp,r1-i,"fwrite");
	}

	if(j<r2){
		tmp=fwrite(a2+j,sizeof(BaseType),r2-j,fout);
		assertIntEQ(tmp,r2-j,"fwrite");
	}

	assertIntEQ(fclose(fin1),0,"fclose");
	assertIntEQ(fclose(fin2),0,"fclose");
	assertIntEQ(fclose(fout),0,"fclose");

	assertIntEQ(unlink(args->fin1),0,args->fin1);
	assertIntEQ(unlink(args->fin2),0,args->fin2);

	pthread_exit(NULL);
}

int partition(BaseType* a, int l, int r) {
   int i=l;
   int j=r+1;
   BaseType t;

   // l sera el pivote
   // y contendra la mediana de l, r y (l+r)/2
   int mid = (l+r)/2;

   if(GT(a[mid],a[r])) {
		 SWAP(a[mid],a[r],t);
   }

   if(GT(a[mid],a[l])) {
		 SWAP(a[mid],a[l],t);
   }

   if(GT(a[l],a[r])) {
		 SWAP(a[l],a[r],t);
	 }

	while (1) {
		do{
			++i;
		}while( !GT(a[i],a[l]) && i <= r );

		do{
			--j;
		}while( GT(a[j],a[l]) && j >= l);

		if( i >= j ) break;

		SWAP(a[i],a[j],t)
	}

	SWAP(a[l],a[j],t)

	return j;
}

int QsortC(BaseType* a, int l,int r) {
   int j;

	if( l < r ) {
 	// divide and conquer
       j = partition( a, l, r);
       //  j=(l+r)/2;
       QsortC( a, l, j-1);
       QsortC( a, j+1, r);
   }
   return 0;
}

void *PQSort(void* a){
	PQSortArgs *args=(PQSortArgs*)a;
	if(args->nth>1){
		int tmp;
		int j = partition(args->a,args->l,args->r);
		int np=1;
		if(args->r - args->l > 0)
			np = (args->nth*(j-args->l))/(args->r-args->l);
		if(np<1) np=1;
		if(np>=args->nth) np=args->nth-1;
		//printf("%d\t%d (%d)\t%d (%d)\n",args->r-args->l,j-args->l,np,args->r-j,args->nth-np);

		pthread_t* th = (pthread_t*)calloc(2,sizeof(pthread_t));
		assertNotNull(th,"calloc");

		PQSortArgs* nargs = (PQSortArgs*)calloc(2,sizeof(PQSortArgs));
		assertNotNull(args,"calloc");

		nargs[0].a=args->a;
		nargs[0].l=args->l;
		nargs[0].r=j-1;
		nargs[0].nth=np;
		tmp=pthread_create(th,NULL,PQSort,(void*)(nargs));
		assertIntEQ(tmp,0,"pthread_create");

		nargs[1].a=args->a;
		nargs[1].l=j+1;
		nargs[1].r=args->r;
		nargs[1].nth=args->nth-np;
		tmp=pthread_create(th+1,NULL,PQSort,(void*)(nargs+1));
		assertIntEQ(tmp,0,"pthread_create");

		tmp=pthread_join(th[0],NULL);
		assertIntEQ(tmp,0,"pthread_join");
		tmp=pthread_join(th[1],NULL);
		assertIntEQ(tmp,0,"pthread_join");

		free(th);
		free(nargs);
	}else{
		QsortC(args->a,args->l,args->r);
	}
	pthread_exit(NULL);
}

unsigned long timestart(){
	struct timeval tv;

	gettimeofday(&tv,NULL);

	return (tv.tv_usec/1000) + (tv.tv_sec*1000);
}

unsigned long timestop(unsigned long start){
	struct timeval tv;

	gettimeofday(&tv,NULL);

	return (tv.tv_usec/1000) + (tv.tv_sec*1000) - start;
}

int psort(int maxsize, int nproc, char* ifile, char* ofile){
	int tmp;
	int max=maxsize;
	int np=nproc;
	if(np<1) np=1;

	printf("Allocating %lu bytes.\n",max*sizeof(BaseType));
	BaseType* a = (BaseType*)calloc(max,sizeof(BaseType));
	assertNotNull(a,"calloc");

	FILE* fin=fopen(ifile,"rt");
	assertNotNull(fin,ifile);

	printf("Stage1: Quicksorts\n");
	unsigned long t = timestart();
	//Read + Quicksort + Write:
	char fname[strlen(ofile)+10];
	int nfile=0;

	do {
		//Read:
		// printf("Reading...\n");
		int n=fread(a,sizeof(BaseType),max,fin);
		if(n==0) break;

		//Quicksort:
		// printf("Quicksort %d\n",n);
		pthread_t th;
		PQSortArgs args;

		args.a=a;
		args.l=0;
		args.r=n-1;
		args.nth=np;
		tmp=pthread_create(&th,NULL,PQSort,(void*)(&args));
		assertIntEQ(tmp,0,"pthread_create");

		//Wait:
		tmp=pthread_join(th,NULL);
		assertIntEQ(tmp,0,"pthread_join");

		//Write:
		// printf("Writing...\n");
		sprintf(fname,"%s.%06d",ofile,nfile);
		FILE* fout=fopen(fname,"wb");
		assertNotNull(fout,"fname");

		tmp=fwrite(a,sizeof(BaseType),n,fout);
		assertIntEQ(tmp,n,"fwrite");

		tmp=fclose(fout);
		assertIntEQ(tmp,0,"fclose");

		nfile++;
	}while(!feof(fin));

	tmp=fclose(fin);
	assertIntEQ(tmp,0,"fclose");

	free(a);

	tmp=unlink(ifile);
	assertIntEQ(tmp,0,"unlink");

	printf("Stage1: %lu\n\n",timestop(t));

	//Merge:
	printf("Stage2: Merge\n");
	t = timestart();

	int curf=0;
	int i=0;
	int j=0;
	pthread_t th[np];
	PMergeArgs args[np];

	curf=0;
	while(nfile-curf>1){
		j=0;
		while(curf<nfile-1 && j<np){
			sprintf(args[j].fin1,"%s.%06d",ofile,curf);
			sprintf(args[j].fin2,"%s.%06d",ofile,curf+1);
			sprintf(args[j].fout,"%s.%06d",ofile,nfile+j);

			// printf("Merge(%d, %d) -> %d\n",curf,curf+1,nfile+j);
			tmp=pthread_create(th+j,NULL,PMerge,(void*)(args+j));
			assertIntEQ(tmp,0,"pthread_create");
			curf+=2;
			j++;
		}

		for(i=0;i<j;i++){
			tmp=pthread_join(th[i],NULL);
			assertIntEQ(tmp,0,"pthread_join");
		}

		nfile+=j;
	}

	printf("Stage2: %lu\n\n",timestop(t));

	//Bin2Text:
	printf("Stage3: Write output file\n");
	t = timestart();

	sprintf(fname,"%s.%06d",ofile,curf);
	if(nfile>0){
		tmp=rename(fname,ofile);
		assertIntEQ(tmp,0,"rename");
	} else {
		FILE *fout=fopen(ofile, "wb");
		assertNotNull(fout,"fopen");
		fclose(fout);
	}

	printf("Stage3: %lu\n",timestop(t));

	return 0;

}
