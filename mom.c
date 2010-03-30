# include <u.h>
# include <libc.h>
# include <bio.h>
<<<<<<< HEAD
# include <stdio.h>
=======
>>>>>>> origin

enum{
	MAXCOL = 1024,	/* maximum number of columns */
	MAXPRES = 15		/* maximum precision */
};
typedef struct Moments Moments;
struct Moments{
	double m1;	/* first moment (mean) */
	double m2;	/* second moment (variance) */
	double m3;	/* third moment (skewness) */
	double m4;	/* fourth moment (kurtosis) */
	long n;		/* number of values seen */
};

Moments momtab[MAXCOL];
int ncol;		/* number of columns */

<<<<<<< HEAD
void fill(char *, Biobuf *);
=======
void fill(Biobuf *);
int iszero(char *);
>>>>>>> origin
void addval(Moments *, double);
void usage(void);

/* use the first 4 moments of a data set to find the mean, standard deviation,
skewness and kurtosis */
void
main(int argc, char *argv[])
{
	Biobuf bstdin, *bf;
	int i, pres;
	char fmt1[7], fmt2[7];
	double sd[MAXCOL];
	long n[MAXCOL];

	pres = 5;
	ARGBEGIN{
	case 'p':
		pres = atoi(EARGF(usage()));
		if(pres < 0 || pres > MAXPRES){
			fprint(2, "%s: invalid precision: %d\n", argv0, pres);
			exits("pres");
		}
		break;
	default:
		usage();
	}ARGEND
	memset(momtab, 0, sizeof(momtab));
<<<<<<< HEAD
	if(argc){
		for(; argc--; argv++)
			if((bf = Bopen(*argv, OREAD)) == nil)
				fprint(2, "%s: unable to open %s: %r\n", argv0, *argv);
			else{
				fill(*argv, bf);
				Bterm(bf);
			}
	}else{
		Binit(&bstdin, 0, OREAD);
		fill("<stdin>", &bstdin);
	}
	sprint(fmt1, "%%.%dg ", pres);
	sprint(fmt2, "%%.%dg\n", pres);
=======
	if(argc)
		while(argc--){
			if((bf = Bopen(*argv++, OREAD)) == nil){
				fprint(2, "%s: unable to open %s: %r\n", argv0, *(argv-1));
				continue;
			}
			fill(bf);
			Bterm(bf);
		}
	else{
		Binit(&bstdin, 0, OREAD);
		fill(&bstdin);
	}
	sprint(fmt1, "%%.%de ", pres);
	sprint(fmt2, "%%.%de\n", pres);
>>>>>>> origin
	for(i = 0; i < ncol; i++)
		print(i < ncol-1 ? fmt1 : fmt2, momtab[i].m1);
	for(i = 0; i < ncol; i++){
		n[i] = momtab[i].n == 1 ? 1 : momtab[i].n-1;
		sd[i] = sqrt(momtab[i].m2 / n[i]);
		print(i < ncol-1 ? fmt1 : fmt2, sd[i]);
	}
	for(i = 0; i < ncol; i++){
<<<<<<< HEAD
		if(sd[i] > 0.0)
=======
		if(sd[i] > 0)
>>>>>>> origin
			print(i < ncol-1 ? fmt1 : fmt2, momtab[i].m3/(n[i]*pow(sd[i], 3)));
		else
			print(i < ncol-1 ? "NaN " : "NaN\n");
	}
	for(i = 0; i < ncol; i++){
<<<<<<< HEAD
		if(sd[i] > 0.0)
=======
		if(sd[i] > 0)
>>>>>>> origin
			print(i < ncol-1 ? fmt1 : fmt2, momtab[i].m4/(n[i]*pow(sd[i], 4)));
		else
			print(i < ncol-1 ? "NaN " : "NaN\n");
	}
	exits(nil);
}

<<<<<<< HEAD
/* fill momtab with data in bf */
void
fill(char *name, Biobuf *bf)
{
	char *line, *tok;
	long nl;
	int ntok;
	double x;

	/* this will fail for lines longer than the bio buffer; Brdstr is an alternative */
	for(nl = 1; line = Brdline(bf, '\n'); nl++){
		line[Blinelen(bf)-1] = '\0';
		ntok = 0;
		while(ntok < MAXCOL && (tok = strtok(ntok == 0 ? line : nil, " \t\r")) != nil){
			if(sscanf(tok, "%lf", &x) == 1)
				addval(momtab+ntok, x);
			ntok++;
		}
		if(ntok > ncol)
			ncol = ntok;
	}
	if(Blinelen(bf)){
		fprint(2, "%s: %s: line %ld too long\n", argv0, name, nl);
		exits("Brdline");
	}
=======
/* fill momtab with data in file bf */
void
fill(Biobuf *bf)
{
	char *tokens[MAXCOL];
	char *line;
	int i, ntok;
	double x;

	/* this will fail for lines longer than the bio buffer; Brstr is an alternative */
	while(line = Brdline(bf, '\n')){
		line[Blinelen(bf)-1] = '\0';
		ntok = 0;
		while(ntok < MAXCOL && (tokens[ntok] = strtok(ntok == 0 ? line : nil, " \t")) != nil)
			ntok++;
		if(ntok > ncol)
			ncol = ntok;
		for(i = 0; i < ntok; i++){
			if((x = atof(tokens[i])) == 0.0 && !iszero(tokens[i]))
				continue;
			addval(momtab+i, x);
		}
	}
}

/* check if s contains a valid zero representation */
int
iszero(char *s)
{
	if(*s == '-' || *s == '+')
		s++;
	return *s == '0';
>>>>>>> origin
}

/* add val to our calculations; formulas deduced from "Formulas for Robust,
One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical
Moments." Philippe P. Pébay, Technical Report SAND2008-6212, Sandia National
Laboratories, September 2008.  */
void
addval(Moments *m, double val)
{
	double n, m1, m2, m3;

	n = m->n;
	m1 = m->m1;
	m2 = m->m2;
	m3 = m->m3;
	m->m1 += (val - m1)/(n+1);
	m->m2 += n*pow(val - m1, 2)/(n+1);
	m->m3 += n*(n-1)*pow(val - m1, 3)/pow(n+1, 2) - 3.0*m2*(val-m1)/(n+1);
	m->m4 += n*(n*n-n+1)*pow(val - m1, 4)/pow(n+1, 3) + 6.0*m2*pow(val - m1, 2)/pow(n+1, 2) - 4.0*m3*(val - m1)/(n+1);
	m->n++;
}

void
usage()
{
<<<<<<< HEAD
	fprint(2, "usage: %s [-p N] [file ...]\n", argv0);
=======
	fprint(2, "usage: %s [-p N] [file] ...\n", argv0);
>>>>>>> origin
	exits("usage");
}
