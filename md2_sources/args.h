#ifndef ARGS_H
#define ARGS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "chem.h"
#include "ion.h"

#define IGNORE_BADWORDS 0
#define TRAP_BADWORDS   1
#define NO_RECURSION    2
#define MAXCHAR		255

void arg_handler (int argc, char * argv[], short bwc);
void parseline (char * ln, char * argv[], int * argc);
void echo_args (int argc, char * argv[]);
void filenames_initialize (void);
void filenames_apply (void);
void filenames_echo (void);
void str_q2n (char *, int);
void str_p2n (char *, int);
#endif
