#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


void parseline (char * ln, char * argv[], int * argc)
{
    char *c;
    int i=-1;
    c=ln;
    *argc=0;
    while (*c)
    {
        if (*c&&!isspace(*c))
            argv[*argc][++i]=*c;
        else 
        {
	    argv[*argc][++i]='\0';
            while(*c&&isspace(*c)) c++;
            if (*c) argv[++*argc][i=0]=*c;
        }
        if (*c) c++;
    }
    ++*argc;
}

char ln[255];
#define MAXCOL 12
#define INIT_STR ("000000000000000000000000000000000000000")
char * buf [MAXCOL]=
{
    INIT_STR, INIT_STR, INIT_STR, 
    INIT_STR, INIT_STR, INIT_STR, 
    INIT_STR, INIT_STR, INIT_STR, 
    INIT_STR, INIT_STR, INIT_STR
};

int main (int argc, char * argv[])
{
    int i=0, n, j;
    short searching;
    double ival, lival, tival;
    FILE * fp = fopen(argv[1], "r");
    i=atoi(argv[2]);
    tival=(double)atof(argv[3]);
    j=atoi(argv[4]);

/*    for (n=0;n<argc;n++) {printf("(%i)%s ", n, argv[n]);}printf("\n");
    printf("file %s i %i tival %.5le j %i\n", argv[1], i, tival, j);
*/
    lival=-999.999;
    searching=1;
    while (searching&&fgets(ln, 255, fp))
    {
	
	if (ln[0]!='#')
	{
	    parseline(ln, buf, &n);
	    ival=atof(buf[i]);
/*	    printf("%f\n", ival); */
	    if (lival!=-999.999)
	    {
		/* we want lival and ival to bracket tival */
		if ((ival<tival&&lival>=tival)||
		    (ival>=tival&&lival<tival))
		{
		    printf("%s\n", buf[j]);
		    searching=0;
		}
	    }
	    lival=ival;
	}
    }

    fclose(fp);    
}
