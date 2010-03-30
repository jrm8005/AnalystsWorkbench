# include <u.h>
# include <stdio.h>
# include <libc.h>
# include <ctype.h>
# include <bio.h>

enum {
	MAXCOLS = 1024,	/* maximum number of columns */
	MAXPRES = 15,		/* maximum output precision */
	INIT = 10,			/* initial allocated size */
	GROW = 2,		/* grow factor */
	IN,
	OUT
};

struct His{
	int dim;			/* histogram dimension */
	int nbins[2];		/* number of bins */
	double *range[2];	/* bin ranges */
	uint *bins;
} his = {
	0,
	{0, 0},
	{nil, nil},
	nil
};

char *status;
char error[] = "Brdline";

char *builds(Biobuf *);
void createhis(char *);
void fill(Biobuf *);
void addvalue(double []);
int binsearch(double, double *, int);
int bincmp(double, double, double);
void printhis(Biobuf *, char *, int, ...);
void *emalloc(ulong);
void usage(void);

void
main(int argc, char *argv[])
{
	Biobuf *bp, bstdin, bstdout;
	char *name, *s, fmt[7];
	int i, j, n, m, pres;
	double *p, *q;
	uint *b;

	SET(s);
	name = nil;
	pres = 5;
	ARGBEGIN{
	case 'f':
		name = EARGF(usage());
		break;
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
	if(name){
		if((bp = Bopen(name, OREAD)) == nil)
			sysfatal("%s: %r", name);
		s = builds(bp);
		Bterm(bp);
	}else if(argc){
		s = *argv++;
		argc--;
	}else
		usage();
	createhis(s);
	if(name)
		free(s);
	if(argc){
		while(argc){
			if((bp = Bopen(*argv, OREAD)) != nil){
				fill(bp);
				Bterm(bp);
			}else
				fprint(2, "%s: %s: %r\n", argv0, *argv);
			argv++;
			argc--;
		}
	}else{
		Binit(&bstdin, 0, OREAD);
		fill(&bstdin);
	}
	Binit(&bstdout, 1, OWRITE);
	sprint(fmt, "%%.%dg ", pres);
	n = his.nbins[0];
	m = his.nbins[1];
	p = his.range[0];
	q = his.range[1];
	b = his.bins;
	if(his.dim == 1)
		for(i = 0; i < n; i++)
			printhis(&bstdout, fmt, 3, p[i], p[i+1], b[i]);
	else
		for(i = 0; i < n; i++)
			for(j = 0; j < m; j++)
				printhis(&bstdout, fmt, 5, p[i], p[i+1], q[j], q[j+1], b[i+j*n]);
	Bterm(&bstdout);
	exits(status);
}

/* builds: build histogram range string */
char *
builds(Biobuf *bp)
{
	char *s, *line;
	int alloc, used, n;

	s = emalloc(INIT*sizeof(char));
	*s = '\0';
	alloc = INIT;
	used = 1;
	for(; line = Brdstr(bp, '\n', 0); free(line)){
		n = Blinelen(bp);
		while(n+1 > alloc-used){
			if((s = realloc(s, GROW*alloc*sizeof(char))) == nil){
				fprint(2, "%s: out of memory\n", argv0);
				exits("nomem");
			}
			alloc *= GROW;
		}
		strcat(s, line);
		used += n;
	}
	return s;
}

/* createhis: build histogram ranges from s */
void
createhis(char *s)
{
	char *t[2], *p;
	double val[2], *r, del;
	int i, j, nbins, nrange, state;

	t[0] = strtok(s, ",");
	t[1] = strtok(nil, ",");
	his.dim = t[1] != nil ? 2 : 1;
	for(i = 0; i < his.dim; i++){
		if(strchr(t[i], ':')){
			if(sscanf(t[i], "%lg:%lg:%d", val, val+1, &nbins) != 3 || nbins < 1){
				fprint(2, "%s: invalid bin range: '%s'\n", argv0, t[i]);
				exits("badpar");
			}
			del = (val[1]-val[0]) / nbins;
			nbins += 2;	/* include under and overflow */
			his.nbins[i] = nbins;
			his.range[i] = emalloc((nbins+1)*sizeof(double));
			r = his.range[i];
			r[0] = Inf(-1);
			for(j = 1; j < nbins; j++)
				r[j] = (j-1)*del + val[0];
			r[nbins] = Inf(1);
		}else{
			state = OUT;
			for(p = t[i], nrange = 0; *p != '\0'; p++)
				if(!isspace(*p) && state == OUT){
					nrange++;
					state = IN;
				}else if(isspace(*p) && state == IN)
					state = OUT;
			if((nbins = nrange-1) < 1){
				fprint(2, "%s: invalid bin range: '%s'\n", argv0, t[i]);
				exits("badpar");
			}
			nbins += 2;	/* include under and overflow */
			his.nbins[i] = nbins;
			his.range[i] = emalloc((nbins+1)*sizeof(double));
			r = his.range[i];
			r[0] = Inf(-1);
			for(j = 1, p = strtok(t[i], " \r\t\n"); j < nbins; j++, p = strtok(nil, " \r\t\n"))
				if(sscanf(p, "%lg", r+j) != 1){
					fprint(2, "%s: invalid bin range: '%s'\n", argv0, t[i]);
					exits("badpar");
				}
			r[nbins] = Inf(1);
		}
	}
	nbins = his.dim == 1 ? his.nbins[0] : his.nbins[0]*his.nbins[1];
	his.bins = emalloc(nbins*sizeof(uint));
	memset(his.bins, 0, nbins*sizeof(uint));
}

/* fill: fill histogram with data in bp */
void
fill(Biobuf *bp)
{
	double val[2];
	char *line;

	while(line = Brdline(bp, '\n')){
		line[Blinelen(bp)-1] = '\0';
		if(sscanf(line, "%lg%lg", val, val+1) == his.dim)
			addvalue(val);
	}
	if(Blinelen(bp))
		status = error;
}

/* addvalue: add val to histogram */
void
addvalue(double val[])
{
	int i, j, n;

	j = 0;
	n = his.nbins[0];
	i = binsearch(val[0], his.range[0], his.nbins[0]);
	if(his.dim == 2)
		j = binsearch(val[1], his.range[1], his.nbins[1]);
	his.bins[n*j + i]++;
}

/* binsearch: find the right bin for val */
int
binsearch(double val, double *v, int n)
{
	int low, high, mid, cond;

	low = 1;			/* skip -Inf */
	high = n-2;		/* and +Inf */
	while(low <= high){
		mid = (low+high) / 2;
		if((cond = bincmp(val, v[mid], v[mid+1])) < 0)
			high = mid-1;
		else if(cond > 0)
			low = mid+1;
		else
			return mid;
	}
	if(val < v[1])
		return 0;
	return n-1;
}

/* bincmp: check if val belongs to min ≤ ξ < max */
int
bincmp(double val, double min, double max)
{
	if(val < min)
		return -1;
	if(val >= max)
		return 1;
	return 0;
}

/* printhis: write histogram to stdout */
void
printhis(Biobuf *bp, char *f, int n, ...)
{
	va_list ap;
	int i;
	double range[4];
	uint bin;
	char str[80];
	char num[25];

	*str = '\0';
	va_start(ap, n);
	for(i = 0; i < n-1; i++){
		range[i] = va_arg(ap, double);
		if(isInf(range[i], -1))
			strcat(str, "underflow ");
		else if(isInf(range[i], 1))
			strcat(str, "overflow ");
		else{
			sprint(num, f, range[i]);
			strcat(str, num);
		}
	}
	bin = va_arg(ap, uint);
	sprint(num, "%d\n", bin);
	strcat(str, num);
	Bprint(bp, str);
	va_end(ap);
}

/* emalloc: allocate memory or die */
void *
emalloc(ulong n)
{
	void *p;

	if((p = malloc(n)) == nil){
		fprint(2, "%s: out of memory\n", argv0);
		exits("nomem");
	}
	return p;
}

void
usage()
{
	fprint(2, "usage: %s [-p N] [-f ranges | 'ranges'] [file ...]\n", argv0);
	exits("usage");
}
