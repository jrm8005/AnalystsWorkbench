# include <u.h>
# include <libc.h>
# include <bio.h>
# include <stdio.h>
# include <ctype.h>

enum {
	MAXPRES = 15,		/* maximum output precision */
	INIT = 3,			/* initial allocated size */
	GROW = 2,		/* grow factor */
};

int dim;			/* histogram dimension */
ulong nentries;		/* total number of entries */
int nbins[3];		/* number of bins for each range */
double *range[3];	/* ranges */
uint *freq;			/* bin frequencies */
int *facs;			/* index factors */
char fmt[6];		/* output format */

int denflag;	/* density flag */
int logflag;	/* logarithm flag */

void creatran(char *, char *, Biobuf *);
void addran(double, char);
void fill(Biobuf *);
void addpoint(double []);
int binsearch(double, double *, int);
int bincmp(double, double, double);
void printhis(Biobuf *, int);
int printr(double, double);
void *emalloc(ulong);
void usage(void);

/* fill a histogram of numeric values */
void
main(int argc, char *argv[])
{
	Biobuf *bp, bstdin, bstdout;
	char *name;
	int i, pres;

	name = nil;
	pres = 5;
	ARGBEGIN{
	case 'd':
		denflag = 1;
		break;
	case 'f':
		name = EARGF(usage());
		break;
	case 'l':
		logflag = 1;
		break;
	case 'p':
		pres = atoi(EARGF(usage()));
		if(pres < 1 || pres > MAXPRES)
			sysfatal("invalid precision: %d", pres);
		break;
	default:
		usage();
	}ARGEND
	if(name){
		if((bp = Bopen(name, OREAD)) == nil)
			sysfatal("%s: %r", name);
		creatran(nil, name, bp);
		Bterm(bp);
	}else if(argc){
		creatran(*argv++, nil, nil);
		argc--;
	}else
		usage();
	if(argc){
		for(i = 0; i < argc; i++, argv++)
			if((bp = Bopen(*argv, OREAD)) != nil){
				fill(bp);
				Bterm(bp);
			}else{
				fprint(2, "%s: %s: %r\n", argv0, *argv);
				if(argc == 1)
					exits("badfile");
			}
	}else{
		Binit(&bstdin, 0, OREAD);
		fill(&bstdin);
	}
	Binit(&bstdout, 1, OWRITE);
	sprint(fmt, "%%.%dg", pres);
	printhis(&bstdout, 0);
	Bterm(&bstdout);

	exits(nil);
}

/* creatran: build histogram ranges from s or bp */
void
creatran(char *s, char *name, Biobuf *bp)
{
	char *p, *ep, *line, *t;
	int i, nl, nx;
	ulong nb;	/* may overflow without too much effort */
	double x, min, max, delt;

	if(s){
		if((t = strdup(s)) == nil)
			sysfatal("out of memory");
		while(p = strtok(dim == 0 ? t : nil, ",")){
			if(sscanf(p, "%lg:%lg:%ld", &min, &max, &nb) != 3 || nb == 0)
				sysfatal("%s: invalid range", p);
			delt = (max-min)/nb;
			for(i = 0; i < nb; i++)
				addran(min + i*delt, 'w');
			addran(max, 'w');	/* remember finite precision */
			addran(0.0, 'l');
		}
		free(t);
	}else
		for(nl = 1; line = Brdstr(bp, '\n', 1); free(line), nl++){
			if(p = strchr(line, '#'))
				*p = '\0';
			for(nx = 0, p = line; *p != '\0'; nx++){
				x = strtod(p, &ep);
				if(ep == p){
					if(!isspace(*p))
						sysfatal("%s:%d: invalid range", name, nl);
					*p = '\0';
					break;
				}
				p = ep;
				addran(x, 'w');
			}
			if(*p == '\0' && p != line){
				if(nx < 2)
					sysfatal("%s:%d: invalid range", name, nl);
				addran(0.0, 'l');
			}
		}
	facs = emalloc(dim*sizeof(int));
	facs[0] = 1;
	for(i = 1; i < dim; i++)
		facs[i] = nbins[i-1]*facs[i-1];
	nb = nbins[0];
	for(i = 1; i < dim; i++)
		nb *= nbins[i];
	if((freq = calloc(nb, sizeof(uint))) == nil)
		sysfatal("out of memory");
}

/* addran: add range element */
void
addran(double x, char c)
{
	static int nalloc[3];
	static double *r;
	int i;

	if(dim > 2)
		sysfatal("too many dimensions");
	/* remember the ±Inf at the extremes */
	if(range[dim] == nil){
		nalloc[dim] = INIT;
		range[dim] = emalloc(nalloc[dim]*sizeof(double));
		*range[dim] = Inf(1);
		r = range[dim] + 1;
		nbins[dim] = 1;
	}
	if(c == 'w'){
		i = r-range[dim];
		if(i >= nalloc[dim]-1){
			range[dim] = realloc(range[dim], GROW*nalloc[dim]*sizeof(double));
			if(range[dim] == nil)
				sysfatal("out of memory");
			r = range[dim] + nalloc[dim]-1;
			nalloc[dim] *= GROW;
		}
		*r++ = x;
		nbins[dim]++;
	}else if(c == 'l'){
		*r = Inf(1);
		dim++;
	}
}

/* fill: fill histogram with data in bp */
void
fill(Biobuf *bp)
{
	static double *pnt = nil;
	double *q;
	char *line, *p;

	if(pnt == nil)
		pnt = emalloc(dim*sizeof(double));
	while(line = Brdline(bp, '\n')){
		line[Blinelen(bp)-1] = '\0';
		if(p = strchr(line, '#'))
			*p = '\0';
		if(p = strtok(line, " \t")){
			q = pnt;
			do{
				if(sscanf(p, "%lg", q) == 1)
					q++;
			}while(p = strtok(nil, " \t"));
			if(q-pnt == dim)
				addpoint(pnt);
		}
	}
}

/* addpoint: add pt to histogram */
void
addpoint(double pnt[])
{
	int i, n;

	n = 0;
	for(i = 0; i < dim; i++)
		n += facs[i]*binsearch(pnt[i], range[i], nbins[i]);
	freq[n]++;
	nentries++;
}

/* binsearch: find the right bin for x */
int
binsearch(double x, double *v, int n)
{
	int low, high, mid, cond;

	low = 1;			/* skip -Inf */
	high = n-2;		/* and +Inf */
	while(low <= high){
		mid = (low+high) / 2;
		if((cond = bincmp(x, v[mid], v[mid+1])) < 0)
			high = mid-1;
		else if(cond > 0)
			low = mid+1;
		else
			return mid;
	}
	if(x < v[1])
		return 0;
	return n-1;
}

/* bincmp: check if min ≤ x < max */
int
bincmp(double x, double min, double max)
{
	if(x < min)
		return -1;
	if(x >= max)
		return 1;
	return 0;
}

char buf[140];
char *bufp = buf;
int ind;
/* printhis: write histogram to stdout */
void
printhis(Biobuf *bp, int n)
{
	int i, k;
	double x;
	char *fm;

	for(i = 0; i < nbins[n]; i++){
		k = printr(*(range[n] + i), *(range[n] + i+1));
		bufp += k;
		ind += i*facs[n];
		if(n < dim-1)
			printhis(bp, n+1);
		else{
			fm = fmt;
			x = freq[ind];
			if(denflag && nentries > 0)
				x /= nentries;
			if(logflag && x == 0.0)
				fm = "-∞";
			else if(logflag)
				x = log10(x);
			sprint(bufp, fm, x);
			Bprint(bp, "%s\n", buf);
		}
		bufp -= k;
		ind -= i*facs[n];
	}
}

/* printr: print range limits */
int
printr(double x, double y)
{
	int n;
	char *fm;

	fm = fmt;
	if(isInf(x, 1))
		fm = "-∞";
	n = sprint(bufp, fm, x);
	*(bufp + n++) = ' ';
	fm = fmt;
	if(isInf(y, 1))
		fm = "+∞";
	n += sprint(bufp+n, fm, y);
	*(bufp + n++) = ' ';
	return n;
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

void
usage()
{
	fprint(2, "usage: %s [-dl] [-p N] [-f ranges | 'ranges'] [file ...]\n", argv0);
	exits("usage");
}
