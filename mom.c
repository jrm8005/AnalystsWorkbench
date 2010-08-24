# include <u.h>
# include <libc.h>
# include <bio.h>
# include <stdio.h>

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

void fill(char *, Biobuf *);
void addval(Moments *, double);
void usage(void);

/* use the first 4 moments of a data set to find the mean, standard deviation,
skewness and kurtosis */
void
main(int argc, char *argv[])
{
	Biobuf bstdin, bstdout, *bf;
	int i, pres;
	char fmt1[7], fmt2[21], *p;
	double sd, skew, kurt;
	long n;

	pres = 5;
	ARGBEGIN{
	case 'p':
		p = EARGF(usage());
		pres = atoi(p);
		if(pres < 1 || pres > MAXPRES)
			sysfatal("invalid precision: %s", p);
		break;
	default:
		usage();
	}ARGEND
	memset(momtab, 0, sizeof(momtab));
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
	sprint(fmt2, "%%.%dg %%.%dg %%.%dg\n", pres, pres, pres);
	Binit(&bstdout, 1, OWRITE);
	for(i = 0; i < ncol; i++){
		Bprint(&bstdout, fmt1, momtab[i].m1);
		if((n = momtab[i].n-1) > 0){
			sd = sqrt(momtab[i].m2 / n);
			skew = momtab[i].m3 / (n*pow(sd, 3));
			kurt = momtab[i].m4 / (n*pow(sd, 4));
			Bprint(&bstdout, fmt2, sd, skew, kurt);
		}else
			Bprint(&bstdout, "NaN NaN NaN\n");
	}
	Bterm(&bstdout);
	exits(nil);
}

/* fill momtab with data in bf */
void
fill(char *name, Biobuf *bf)
{
	char *line, *tok;
	long nl;
	int ntok;
	double x;

	/* this will fail for lines longer than the bio buffer; Brdstr is an alternative */
	for(nl = 1; (line = Brdline(bf, '\n')); nl++){
		line[Blinelen(bf)-1] = '\0';
		ntok = 0;
		while(ntok < MAXCOL && (tok = strtok(ntok == 0 ? line : nil, " \t\r")) ){
			if(sscanf(tok, "%lg", &x) == 1)
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
	fprint(2, "usage: %s [-p N] [file ...]\n", argv0);
	exits("usage");
}
