#include "decs.h"
#include <ctype.h>

// ppm image information:
#define NCOLORS	(256)

/* ppm image variables */
static unsigned int color_map[3][NCOLORS];
static int ppm_height, ppm_ncolors, ppm_width;

int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

#define RED		(0)
#define GREEN		(1)
#define BLUE		(2)
#define TOP_FRAC	(0.005)	/* Brightest and dimmest pixels cut out of color scale */
#define BOT_FRAC	(0.005)

void make_ppm(double p[NX][NY], double freq, char filename[])
{

	int i, j, k, icolor, npixels;
	double *q, min, max, scale, iq, liq;
	FILE *fp;
	static int firstc = 1 ;
	double *alloc_double_array( int ndata )  ;
	void get_color_map(void) ;

	if(firstc) {
		firstc = 0 ;
		get_color_map() ;
	}

	npixels = NX * NY;

	q = alloc_double_array(npixels);
	k = 0 ;
        for (i = 0; i < NX; i++)
	for (j = 0; j < NY; j++) {
		q[k] = p[i][j];
		k++ ;
	}

	qsort(q, npixels, sizeof(double), compare_doubles);
	if (q[0] < 0) {		/* must be log scaling */
		min = q[(int) (npixels * BOT_FRAC)];
	} else {
		i = 0;
		while (i < npixels - 1 && q[i] == 0.)
			i++;
		min = q[i];
	}
	max = q[npixels - (int) (npixels * TOP_FRAC)];

	/*
	min = -0.1 ;
	max = 1.1 ;
	*/

	scale = ((double) ppm_ncolors) / (max - min);

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Failed to open ppm file %s\n", filename);
		exit(125);
	}

	/* write out header information */
	fprintf(fp,
		"P6\n#  min=%g  , max=%g \n#  frequency=%g \n%d %d\n%d\n",
		min, max, freq, NX, NY, ppm_ncolors);
	fflush(fp);

	for (j = NY-1; j >= 0; j--) 
	for (i = 0; i < NX; i++) 
	{
		iq = p[i][j];
		liq = scale * (iq - min);
		icolor = (int) liq;
		if (icolor > ppm_ncolors) {
			icolor = ppm_ncolors;
		}
		if (icolor < 0) {
			icolor = 0;
		}
		fputc(color_map[RED][icolor], fp);
		fputc(color_map[GREEN][icolor], fp);
		fputc(color_map[BLUE][icolor], fp);
	}

	fclose(fp);

	return;
}

#undef TOP_FRAC
#undef BOT_FRAC

/*******************************************************************************/
/*******************************************************************************
   get_color_map():
   ---------------

     -- reads in a plain-type ppm file to generate the color map 
     -- see  http://netpbm.sourceforge.net/doc/ppm.html  or "man ppm" 
          for PPM format specifications 

*******************************************************************************/
void get_color_map(void)
{

	int i, c0, c1, c2, iw, ibody;
	char c, word[100];
	int ipixel, icolor;
	FILE *fp;
	const char magic_raw[] = "P6";
	const char magic_plain[] = "P3";
	enum GET_type { get_magic, get_width, get_height, get_ncolors,
		    get_pixels, get_newline };
	enum GET_type read_type, old_read_type;
	enum PPM_type { ppm_plain, ppm_raw };
	enum PPM_type ppm_type, ppm_body;

	fprintf(stdout, "Starting get_color_map().... \n");
	fflush(stdout);

	if ((fp = fopen("map.ppm", "r")) == NULL) {
		fflush(stderr);
		fprintf(stderr,
			"get_color_map(): Cannot open map.ppm !! \n");
		fflush(stderr);
		exit(1);
	}

	old_read_type = read_type = get_magic;	/* Need to read in the magic number or ppm type first */
	ppm_type = ppm_plain;	/* Always read in the header as plain text */
	iw = 0;
	icolor = ipixel = 0;

	/* Read in the entire stream as a character stream */
	while ((c = getc(fp)) != EOF) {

		//fprintf(stderr,"%c",c);fflush(stderr);

		switch (read_type) {

		case get_magic:
			word[iw++] = c;
			if (iw == 2) {
				word[iw] = '\0';
				if (strcmp(word, magic_raw) == 0) {
					ppm_body = ppm_raw;
				} else if (strcmp(word, magic_plain) == 0) {
					ppm_body = ppm_plain;
				} else {
					fprintf(stderr,
						"get_color_map(): error getting magic number, type = %s \n",
						word);
					fflush(stderr);
					fclose(fp);
					exit(1);
				}
				old_read_type = read_type;
				read_type = get_width;
				iw = 0;
			}
			break;

		case get_width:
			if (c == '#') {
				old_read_type = get_width;
				read_type = get_newline;
				break;
			}
			if (isdigit(c)) {
				word[iw++] = c;
			} else if (iw > 0) {
				word[iw] = '\0';
				ppm_width = atoi(word);
				old_read_type = read_type;
				read_type = get_height;
				iw = 0;
			}
			break;

		case get_height:
			if (c == '#') {
				old_read_type = get_height;
				read_type = get_newline;
				break;
			}
			if (isdigit(c)) {
				word[iw++] = c;
			} else if (iw > 0) {
				word[iw] = '\0';
				ppm_height = atoi(word);
				old_read_type = read_type;
				read_type = get_ncolors;
				iw = 0;
			}
			break;

		case get_ncolors:
			if (c == '#') {
				old_read_type = get_ncolors;
				read_type = get_newline;
				break;
			}
			if (isdigit(c)) {
				word[iw++] = c;
			} else if (iw > 0) {
				word[iw] = '\0';
				ppm_ncolors = atoi(word);
				old_read_type = read_type;
				read_type = get_pixels;
				iw = 0;
				if (ppm_ncolors > NCOLORS) {
					fprintf(stderr,
						"get_color_map(): too many colors: ncolors = %d,  maxcolors = %d \n",
						ppm_ncolors, NCOLORS);
					fflush(stderr);
					fclose(fp);
					exit(1);
				}
			}
			break;

		case get_pixels:
			/* We assume that pixels come in threes and there are no comments within a pixel */
			if (c == '#') {
				old_read_type = get_pixels;
				read_type = get_newline;
				break;
			}
			if (isdigit(c)) {
				word[iw++] = c;
			} else if (iw > 0) {
				word[iw] = '\0';
				color_map[icolor++][ipixel] = atoi(word);
				iw = 0;
				if (icolor == 3) {
					icolor = 0;
					ipixel++;
				}
			}
			break;

		case get_newline:
			if (c == '\n') {
				read_type = old_read_type;
			}
			break;

		}		/* switch off */

	}			/* end while */

	if (ipixel != (ppm_width * ppm_height)) {
		fprintf(stderr,
			"get_color_map(): missing pixels: npixel = %d,  width = %d,  height = %d \n",
			ipixel, ppm_width, ppm_height);
		fflush(stderr);
		fclose(fp);
		exit(1);
	}

	if ((ppm_ncolors + 1) != ipixel) {
		fprintf(stderr,
			"get_color_map(): missing colors: npixel = %d,  ppm_ncolors = %d \n",
			ipixel, ppm_ncolors);
		fflush(stderr);
		fclose(fp);
		exit(1);
	}

	fclose(fp);


  /********************************************************************************
    test get_color_map():  Make sure that testmap.ppm looks identical to your map.ppm file. 
  ********************************************************************************/
	if ((fp = fopen("testmap.ppm", "w")) == NULL) {
		exit(1);
	}
	fprintf(fp, "P6\n%d %d\n%d\n", ppm_width, ppm_height, ppm_ncolors);
	fflush(fp);
	for (i = 0; i <= ppm_ncolors; i++) {
		fputc(color_map[RED][i], fp);
		fputc(color_map[GREEN][i], fp);
		fputc(color_map[BLUE][i], fp);
	}
	fclose(fp);


	fprintf(stdout, "Ending get_color_map().... \n");
	fflush(stdout);

	return;
}

/*********************************************************************************************
  alloc_double_array():
     -- returns with a pointer to an array of size "ndata"
 *********************************************************************************************/
double *alloc_double_array( int ndata ) 
{ 
  double *pa;

  pa = (double *) calloc(ndata,sizeof(double));
  if( pa == NULL ) { 
    fprintf(stderr,"Error allocating a double array of length = %d \n", ndata);
    fprintf(stderr,"....Exiting...\n");
    fflush(stderr);
    exit(1);
  }
  return(pa);
}


#undef RED
#undef GREEN
#undef BLUE
