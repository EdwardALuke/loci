/*
 * getargs.h - include file for getargs.c
 */

#ifndef _GETARGS_H_
#define _GETARGS_H_

#ifdef __cplusplus
extern "C" {
#endif

extern int argind;      /* index into argv array */
extern int argerr;      /* flag for printing errors */
extern char *argarg;    /* pointer to argument string */

void print_usage (
    const char **usgmsg,      /* usage message */
    const char *errmsg        /* error message */
);

int getargs (
    int argc,           /* number of arguments */
    char **argv,        /* argument list */
    char *ostr          /* option list */
);

#ifdef __cplusplus
}
#endif

#endif  /* _GETARGS_H_ */
