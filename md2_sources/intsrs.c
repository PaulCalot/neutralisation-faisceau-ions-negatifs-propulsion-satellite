#include <stdio.h>
#include <string.h>
#include <stdlib.h>
void main (int argc, char * argv[])
{
    int i = 0, j=0;
    int n, f, d;
    short z=1;
    for (i=1;i<argc;i++)
    {
	if (!strcmp(argv[i], "-n")) n = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-f")) f = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-d")) d = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-nz")) z=0;
	else exit(0);
    }
    for (i=0;i<n;i++) 
    {
	j=f+d*i;
	printf("%s%s%s%i ", 
	    (z&&j<1000?"0":""), (z&&j<100?"0":""), (z&&j<10?"0":""), j);
    }
    printf("\n");
}
