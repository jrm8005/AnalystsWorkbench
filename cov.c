# include <u.h>
# include <libc.h>
# include <bio.h>
# include <stdio.h>
# include <ctype.h>

enum{
	MAXPRES = 15		/* maximum precision */
};
typedef struct Csum Csum;
struct Csum{
	double xm1;	/* means */
	double ym1;
	double xm2;	/* variances */
	double ym2;
	double sum;	/* corrected sum */
	long n;		/* number of values seen */
};

Csum **csumtab;
int nfields;
int pearson;	/* Pearson's correlation coefficient flag */

void fill(Biobuf *);
void add(double *);
double get(int, int);
void *emalloc(ulong);
void *ecalloc(ulong, ulong);
void usage(void);

/* find the covariance matrix for a data set */
void
main(int argc, char *argv[])
{
	Biobuf *bp, bstdin, bstdout;
	int i, j, pres;
	char fmt[6];
	char *fm;
	double val;

	pres = 5;
	ARGBEGIN{
	case 'p':
		pres = atoi(EARGF(usage()));
		if(pres < 1 || pres > MAXPRES)
			sysfatal("invalid precision: %d", pres);
		break;
	case 'r':
		pearson = 1;
		break;
	default:
		usage();
	}ARGEND
	if(argc){
		for(; argc--; argv++)
			if((bp = Bopen(*argv, OREAD)) == nil)
				fprint(2, "%s: unable to open %s: %r\n", argv0, *argv);
			else{
				fill(bp);
				Bterm(bp);
			}
	}else{
		Binit(&bstdin, 0, OREAD);
		fill(&bstdin);
	}
	sprint(fmt, "%%.%dg", pres);
	Binit(&bstdout, 1, OWRITE);
	for(i = 0; i < nfields; i++)
		for(j = 0; j < nfields; j++){
			val = get(i, j);
			fm = fmt;
			if(isInf(val, 1))
				fm = "+∞";
			else if(isInf(val, -1))
				fm = "-∞";
			Bprint(&bstdout, fm, val);
			Bprint(&bstdout, j < nfields-1 ? " " : "\n");
		}
	Bterm(&bstdout);
	exits(nil);
}

/* fill: read data in bp */
void
fill(Biobuf *bp)
{
	static double *v = nil;
	double *q;
	char *line, *p, *s;
	Csum **c;
	int n;

	while(line = Brdline(bp, '\n')){
		line[Blinelen(bp)-1] = '\0';
		if(p = strchr(line, '#'))
			*p = '\0';
		if(csumtab == nil){
			if((s = strdup(line)) == nil)
				sysfatal("out of memory");
			for(n = 0; strtok(n == 0 ? s : nil, " \t"); n++)
				;
			if(n){
				nfields = n;
				v = emalloc(nfields*sizeof(double));
				csumtab = emalloc(nfields*sizeof(Csum *));
				for(c = csumtab; n > 0; c++)
					*c = ecalloc(n--, sizeof(Csum));
			}
			free(s);
		}
		if(p = strtok(line, " \t")){
			q = v;
			do{
				if(sscanf(p, "%lg", q) != 1)
					*q = Inf(1);
				if(++q - v >= nfields)
					break;
			}while(p = strtok(nil, " \t"));
			if(q - v == nfields)
				add(v);
		}
	}
}

/* add: add values in v to the calculations; formulas taken from B. P. Welford
(1962)."Note on a method for calculating corrected sums of squares and
products".  Technometrics 4(3):419–420 */
void
add(double *v)
{
	int i, j;
	double x, y, xm1, ym1, xm2, ym2;
	long n;
	Csum *c;

	for(i = 0; i < nfields; i++)
		for(j = i; j < nfields; j++)
			if(!isInf(v[i], 1) && !isInf(v[j], 1)){
				c = *(csumtab+i) + j-i;
				x = v[i];
				y = v[j];
				xm1 = c->xm1;
				ym1 = c->ym1;
				xm2 = c->xm2;
				ym2 = c->ym2;
				n = c->n;
				c->xm1 += (x - xm1)/(n+1);
				c->ym1 += (y - ym1)/(n+1);
				c->xm2 += n*pow(x - xm1, 2)/(n+1);
				c->ym2 += n*pow(y - ym1, 2)/(n+1);
				c->sum += n*(x - xm1)*(y - ym1)/(n+1);
				c->n++;
			}
}

/* get: return value of cov_ij */
double
get(int i, int j)
{
	int aux;
	double xsd, ysd, cov;
	Csum *c;

	if(j < i){	/* it's a symmetric matrix */
		aux = i;
		i = j;
		j = aux;
	}
	c = *(csumtab+i) + j-i;
	if(c->n){
		xsd = sqrt(c->xm2/c->n);
		ysd = sqrt(c->ym2/c->n);
		cov = c->sum/c->n;
	}else
		return Inf(1);
	if(pearson){
		if(xsd > 0.0 && ysd > 0.0)
			cov /= xsd*ysd;
		else if(cov >= 0.0)
			cov = Inf(1);
		else if(cov < 0.0)
			cov = Inf(-1);
	}
	return cov;
}

/* emalloc: allocate memory or die */
void *
emalloc(ulong n)
{
	void *p;

	if((p = malloc(n)) == nil)
		sysfatal("out of memory");
	return p;
}

/* ecalloc: allocate memory or die */
void *
ecalloc(ulong n, ulong size)
{
	void *p;

	if((p = calloc(n, size)) == nil)
		sysfatal("out of memory");
	return p;
}

void
usage()
{
	fprint(2, "usage: %s [-r] [-p N] [file ...]\n", argv0);
	exits("usage");
}
