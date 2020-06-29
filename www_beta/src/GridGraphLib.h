/* 
 * Copyright (C) 1999 Anthony Lomax <lomax@faille.unice.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/*  GridGraphicsLib.h

	include file for grid graphics functions

*/



/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */



/* include files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*#include <errno.h> */
#include <time.h>
/*#include <float.h> */


/* defines */

#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif





/*------------------------------------------------------------/ */
/* macros */
/*------------------------------------------------------------/ */


/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* structures */
/*------------------------------------------------------------/ */


	/* colors */

typedef struct {
    double r;	/* red */
    double g;	/* red */
    double b;	/* red */
} ColRGB;


	/* stations */

struct GrStation {
    int n_sta;      	/* station id number */
    char name[11];      /* name */
    int txt_pos;      	/* text position code 0=lower left, 1=ul, 2=ur, 3=lr */
    double txt_size;   	/* text size as fraction of default size */
    char col[11];		/* color */
    double zloc;    	/* z location */
    double yloc;    	/* y location */
    double xloc;    	/* x location */
};
EXTERN_TXT int gr_num_sta;		/* number of stations read */

#define MAX_NUM_STA 100
EXTERN_TXT struct GrStation sta_array[MAX_NUM_STA];


	/* graphics text */

struct GrTxt {
    char plot[5];		/* view SECT or MAP */
    char text[51];      /* name */
    char col[11];		/* color */
    double size;    	/* text size factor */
    double zloc;    	/* z location */
    double yloc;    	/* y location */
    double xloc;    	/* x location */
    int rot;   		 	/* rotation in deg */
};

EXTERN_TXT int gr_num_txt;		/* number of text strings read */

#define MAX_NUM_TXT 100
EXTERN_TXT struct GrTxt txt_array[MAX_NUM_TXT];

	/* graphics lines */

struct GrLine {
    char plot[5];		/* view SECT or MAP */
    double z1;    		/* z location */
    double y1;    		/* y location */
    double x1;    		/* x location */
    double z2;    		/* z location */
    double y2;    		/* y location */
    double x2;    		/* x location */
    char col[11];		/* color */
    char style[11];     /* style */
};

EXTERN_TXT int gr_num_line;		/* number of lines read */

#define MAX_NUM_LINE 100
EXTERN_TXT struct GrLine line_array[MAX_NUM_LINE];

	/* map line files */

struct GrMapLines {
    char format[MAXLINE];      	/* file format */
    char name[2*MAXLINE];		/* filename */
    char line_style[11];	/* linesytle */
    ColRGB rgb;			/* line color */
};

EXTERN_TXT int gr_num_map_files;	/* number of map files */

#define MAX_NUM_MAP 10
EXTERN_TXT struct GrMapLines mapfile[MAX_NUM_MAP];







/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* globals  */
/*------------------------------------------------------------/ */



/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

void SetGraphicConstants(void);
int ReadGraphicsInput(FILE* , char* , char* , int );
int get_maplines(char* );
char *get_timestamp();
int Err2GMT(FILE* , double , double , double , double , double , double );
Vect2D* Cov2Ellipse(double , double, double, double , double , Vect2D* , int );
int Vect2DArray2GMT(FILE* , Vect2D* , int );
Vect3D* toEllipsoid3D(Vect3D , Vect3D , Vect3D , int );
Vect2D*  Vect3D2To2D(Vect3D *, Vect2D* , int , int );

/* */
/*------------------------------------------------------------/ */


