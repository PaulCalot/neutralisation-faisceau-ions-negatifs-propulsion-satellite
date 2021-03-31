#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

double variableOne;
int variableTwo;

typedef struct VAR
{
    char * name;
    void * val;
    enum {INT, DBL} dt;
} key;

#define IVAR(x)  {#x, &x, INT}
#define DVAR(x)  {#x, &x, DBL}

key keyList[] = 
{
    IVAR(variableTwo), 
    DVAR(variableOne), 
};

void main (int argc, char * argv[])
{
    int i;
    
   *((int*)(keyList[0].val))=atoi(argv[1]);
    *((double*)(keyList[1].val))=atof(argv[2]);
  
    for (i=0;i<2;i++)
    {
	printf("%s ", keyList[i].name);
	if (keyList[i].dt==INT) 
	    printf("%i\n", *((int*)(keyList[i].val)));
	else if (keyList[i].dt==DBL) 
	    printf("%7.4f\n", *((double*)(keyList[i].val)));
    }
}
