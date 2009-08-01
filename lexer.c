#include "lexer.h"

#define chatty 0

/* Routine to interpret input and return estimate of type.
 * This needs to read at most 2 lines
 */
input_type lex()
{
	int i=0, ndigits=0, nlines=0;
	char *c;
	
	while (c = getc(stdin)) {
		if (c == '\n') {
			ndigits=0; 
			++nlines;
			if ( c = getc(stdin) != EOF && nlines < 2) {
				ungetc(c, stdin);
				i = 0;
			}
		}
		
		if (isdigit(c) || c == '.') {
			if (ndigits == 0)
				++i;
			++ndigits;
		}

		
		if (isspace(c)) {
			ndigits=0;
		}
		
		if (nlines == 2 || c == EOF)
			break;
	}
	
	if (chatty) {
		printf("nlines: %d, i: %d\n", nlines, i);
	}
	
	if (nlines == 1) {
		if (i < 2)
			return SCALAR;
		else return VECTOR;
	}
	
	if (i > 1)
		return MATRIX;
	else return VECTOR;
}

/* TODO: Code for detecting special case 2-column matrix. To be integrated later. */
#if 0 /* I know rob would kill me for this */
	while (c = getc(stdin)) {
		if (nlines >= 2 || c == EOF)
			break;
		
		if (c == '\n') {
			printf("\n");
			ndigits=0; 
			++nlines;
			if ( c = getc(stdin) != EOF && nlines < 2) {
				ungetc(c, stdin);
				i = 1;
			}
		}
		
		if (isdigit(c) || c == '.') {
			++ndigits;
		}

		
		if (isspace(c)) {
			ndigits=0;
			++i;
		}	
		printf("i is: %d\n", i);

	}
	
	if (chatty) {
		printf("nlines: %d, i: %d\n", nlines, i);
	}
	
	if (nlines == 1) {
		if (i < 2)
			printf("scalar\n");
		else printf("vector\n");
		return;
	}
	
	if (i >= 1) {
		if (i ==2) {
			printf("2 column matrix\n");
		} else printf("matrix\n");
	}
	else printf("vector\n");
#endif