/*
 * Programmation Parall�le - Avril 2012
 * Polytech'Paris 
 * Universit� Pierre et Marie Curie
 * Calcul de l'ensemble de Mandelbrot, Version s�quentielle
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>	/* chronometrage */
#include <string.h>     /* pour memset */
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#include "rasterfile.h"

char info[] = "\
Usage:\n\
      mandel dimx dimy xmin ymin xmax ymax prof\n\
\n\
      dimx,dimy : dimensions de l'image a generer\n\
      xmin,ymin,xmax,ymax : domaine a calculer dans le plan complexe\n\
      prof : nombre maximale d'iteration\n\
\n\
Quelques exemples d'execution\n\
      mandel 800 800 0.35 0.355 0.353 0.358 200\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 500\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 300\n\
      mandel 800 800 -1.48478 0.00006 -1.48440 0.00044 100\n\
      mandel 800 800 -1.5 -0.1 -1.3 0.1 10000\n\
";

double my_gettimeofday(){
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

/**
 * Convertion entier (4 octets) LINUX en un entier SUN
 * @param i entier � convertir
 * @return entier converti
 */

int swap(int i) {
	int init = i; 
	int conv;
	unsigned char *o, *d;
	  
	o = ( (unsigned char *) &init) + 3; 
	d = (unsigned char *) &conv;

	*d++ = *o--;
	*d++ = *o--;
	*d++ = *o--;
	*d++ = *o--;

	return conv;
}

/*** 
 * Par Francois-Xavier MOREL (M2 SAR, oct2009): 
 */

unsigned char power_composante(int i, int p) {
	unsigned char o;
	double iD=(double) i;

	iD/=255.0;
	iD=pow(iD,p);
	iD*=255;
	o=(unsigned char) iD;
	return o;
}

unsigned char cos_composante(int i, double freq) {
	unsigned char o;
	double iD=(double) i;
	iD=cos(iD/255.0*2*M_PI*freq);
	iD+=1;
	iD*=128;

	o=(unsigned char) iD;
	return o;
}

/*** 
 * Choix du coloriage : definir une (et une seule) des constantes
 * ci-dessous :  
 */
//#define ORIGINAL_COLOR
#define COS_COLOR 

#ifdef ORIGINAL_COLOR
#define COMPOSANTE_ROUGE(i)    ((i)/2)
#define COMPOSANTE_VERT(i)     ((i)%190)
#define COMPOSANTE_BLEU(i)     (((i)%120) * 2)
#endif /* #ifdef ORIGINAL_COLOR */
#ifdef COS_COLOR
#define COMPOSANTE_ROUGE(i)    cos_composante(i,13.0)
#define COMPOSANTE_VERT(i)     cos_composante(i,5.0)
#define COMPOSANTE_BLEU(i)     cos_composante(i+10,7.0)
#endif /* #ifdef COS_COLOR */

/**
 *  Sauvegarde le tableau de donn�es au format rasterfile
 *  8 bits avec une palette de 256 niveaux de gris du blanc (valeur 0)
 *  vers le noir (255)
 *    @param nom Nom de l'image
 *    @param largeur largeur de l'image
 *    @param hauteur hauteur de l'image
 *    @param p pointeur vers tampon contenant l'image
 */

void sauver_rasterfile( char *nom, int largeur, int hauteur, unsigned char *p) {
	FILE *fd;
	struct rasterfile file;
	int i;
	unsigned char o;

	if ( (fd=fopen(nom, "w")) == NULL ) {
		printf("erreur dans la creation du fichier %s \n",nom);
		exit(1);
	}

	file.ras_magic  = swap(RAS_MAGIC);	
	file.ras_width  = swap(largeur);	  /* largeur en pixels de l'image */
	file.ras_height = swap(hauteur);         /* hauteur en pixels de l'image */
	file.ras_depth  = swap(8);	          /* profondeur de chaque pixel (1, 8 ou 24 )   */
	file.ras_length = swap(largeur*hauteur); /* taille de l'image en nb de bytes		*/
	file.ras_type    = swap(RT_STANDARD);	  /* type de fichier */
	file.ras_maptype = swap(RMT_EQUAL_RGB);
	file.ras_maplength = swap(256*3);

	fwrite(&file, sizeof(struct rasterfile), 1, fd); 

	/* Palette de couleurs : composante rouge */
	i = 256;
	while( i--) {
	o = COMPOSANTE_ROUGE(i);
	fwrite( &o, sizeof(unsigned char), 1, fd);
	}

	/* Palette de couleurs : composante verte */
	i = 256;
	while( i--) {
	o = COMPOSANTE_VERT(i);
	fwrite( &o, sizeof(unsigned char), 1, fd);
	}

	/* Palette de couleurs : composante bleu */
	i = 256;
	while( i--) {
	o = COMPOSANTE_BLEU(i);
	fwrite( &o, sizeof(unsigned char), 1, fd);
	}

	// pour verifier l'ordre des lignes dans l'image : 
	//fwrite( p, largeur*hauteur/3, sizeof(unsigned char), fd);

	// pour voir la couleur du '0' :
	// memset (p, 0, largeur*hauteur);

	fwrite( p, largeur*hauteur, sizeof(unsigned char), fd);
	fclose( fd);
}

/**
 * �tant donn�e les coordonn�es d'un point \f$c=a+ib\f$ dans le plan
 * complexe, la fonction retourne la couleur correspondante estimant
 * � quelle distance de l'ensemble de mandelbrot le point est.
 * Soit la suite complexe d�fini par:
 * \f[
 * \left\{\begin{array}{l}
 * z_0 = 0 \\
 * z_{n+1} = z_n^2 - c
 * \end{array}\right.
 * \f]
 * le nombre d'it�rations que la suite met pour diverger est le
 * nombre \f$ n \f$ pour lequel \f$ |z_n| > 2 \f$. 
 * Ce nombre est ramen� � une valeur entre 0 et 255 correspond ainsi a 
 * une couleur dans la palette des couleurs.
 */

unsigned char xy2color(double a, double b, int prof) {
	double x, y, temp, x2, y2;
	int i;

	x = y = 0.;
	for( i=0; i<prof; i++) {
		/* garder la valeur pr�c�dente de x qui va etre ecrase */
		temp = x;
		/* nouvelles valeurs de x et y */
		x2 = x*x;
		y2 = y*y;
		x = x2 - y2 + a;
		y = 2*temp*y + b;
		if( x2 + y2 >= 4.0) break;
	}
	return (i==prof)?255:(int)((i%255)); 
}

/* 
 * Partie principale: en chaque point de la grille, appliquer xy2color
 */

int main(int argc, char *argv[]) {
	/* Domaine de calcul dans le plan complexe */
	double xmin, ymin;
	double xmax, ymax;
	/* Dimension de l'image */
	int w,h;
	/* Pas d'incrementation */
	double xinc, yinc;
	/* Profondeur d'iteration */
	int prof;
	/* Image resultat */
	unsigned char	*ima, *pima;
	/* Variables intermediaires */
	int  i, j;
	double x, y;
	/* Chronometrage */
	double debut, fin;

	int nb_proc;

	/* debut du chronometrage */
	debut = my_gettimeofday();


	if( argc == 1) fprintf( stderr, "%s\n", info);

	/* Valeurs par defaut de la fractale */
	xmin = -2; ymin = -2;
	xmax =  2; ymax =  2;
	w = h = 800;
	prof = 10000;
	nb_proc = 1;

	/* Recuperation des parametres */
	if( argc > 1) nb_proc    = atoi(argv[1]);
	if( argc > 2) w    = atoi(argv[2]);
	if( argc > 3) h    = atoi(argv[3]);
	if( argc > 4) xmin = atof(argv[4]);
	if( argc > 5) ymin = atof(argv[5]);
	if( argc > 6) xmax = atof(argv[6]);
	if( argc > 7) ymax = atof(argv[7]);
	if( argc > 8) prof = atoi(argv[8]);

	/* Calcul des pas d'incrementation */
	xinc = (xmax - xmin) / (w-1);
	yinc = (ymax - ymin) / (h-1);

	/* affichage parametres pour verificatrion */
	fprintf( stderr, "Domaine: {[%lg,%lg]x[%lg,%lg]}\n", xmin, ymin, xmax, ymax);
	fprintf( stderr, "Increment : %lg %lg\n", xinc, yinc);
	fprintf( stderr, "Prof: %d\n",  prof);
	fprintf( stderr, "Dim image: %dx%d\n", w, h);

	/* Allocation memoire du tableau resultat */  
	pima = ima = (unsigned char *)malloc( w*h*sizeof(unsigned char));

	if( ima == NULL) {
		fprintf( stderr, "Erreur allocation m�moire du tableau \n");
		return 0;
	}

	omp_set_num_threads(nb_proc);

	#pragma omp parallel {
		/* Traitement de la grille point par point */
		y = ymin; 
		#pragma omp for private(j,pima,x,y) schedule (runtime)
		for (i = 0; i < h; i++) {
			y = ymin + i * yinc;
			pima = &ima[i*w];	
			x = xmin;
			for (j = 0; j < w; j++) {
				// printf("%d\n", xy2color( x, y, prof));
				// printf("(x,y)=(%g;%g)\t (i,j)=(%d,%d)\n", x, y, i, j);
				*pima++ = xy2color( x, y, prof); 
				x += xinc;
			}
			y += yinc; 
		}
	}

	/* fin du chronometrage */
	fin = my_gettimeofday();
	fprintf( stderr, "Temps total de calcul : %g sec\n", 
	fin - debut);
	fprintf( stdout, "%g\n", fin - debut);

	/* Sauvegarde de la grille dans le fichier resultat "mandel.ras" */
	sauver_rasterfile( "mandel.ras", w, h, ima);

	return 0;
}