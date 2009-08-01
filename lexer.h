#ifndef LEXER_H
#define LEXER_H

#include <stdio.h>

/* Type of input. TAGGED* is for input with type explicitly put on line 0 */
typedef enum input_type input_type;
enum {
	MATRIX,
	VECTOR,
	SCALAR,
	TAGGED_MATRIX,
	TAGGED_VECTOR,
	TAGGED_SCALAR
} input_type;

input_type lex(FILE *f);

#endif LEXER_H
