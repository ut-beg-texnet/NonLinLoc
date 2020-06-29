/*
 * Copyright (C) 1999-2010 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


/* velmod.h include file */

/* 	History:

        see wls.c

 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define MAXFILE 200

/*#define REAL_MIN min_normal()*/
#ifndef REAL_MIN
#define REAL_MIN ((double) 1.0e-30)
#endif
#ifndef REAL_MAX
#define REAL_MAX ((double) 1.0e+30)
#endif
#ifndef REAL_EPSILON
#define REAL_EPSILON ((double) 1.0e-6)
#endif

#ifndef LARGE_FLOAT
#define LARGE_FLOAT 1.0e20F
#endif

/* the following defines needed for some old versions of cc */
#ifndef SEEK_SET
#define SEEK_SET 0
#endif
#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif
#ifndef REAL_MAX
#define REAL_MAX max_normal()
#endif

#ifndef MAXLINE
#define MAXLINE 101
#endif

#ifndef MAXLINE_LONG
#define MAXLINE_LONG 4*MAXLINE
#endif

/* externally defined names */
extern int prog_mode_3d; /* 0 = 2D, 1 = 3D calculation */
//extern void puterr(char* );

EXTERN_TXT double min_x_cut; /* minimum x distance cutoff */
EXTERN_TXT double max_x_cut; /* maximum x distance cutoff */
EXTERN_TXT double min_y_cut; /* minimum y distance cutoff */
EXTERN_TXT double max_y_cut; /* maximum y distance cutoff */
EXTERN_TXT double min_z_cut; /* minimum y distance cutoff */
EXTERN_TXT double max_z_cut; /* maximum y distance cutoff */


/**************************************************************************/
/* model definitions */
/**************************************************************************/


/* 2D to 3D transformation */

EXTERN_TXT int prog_mode_Mod2D3D; /* 1 = 3D calculation with 2D model */
EXTERN_TXT double Mod2D3D_origx, Mod2D3D_origy, Mod2D3D_rot;
EXTERN_TXT double Mod2D3D_cosang, Mod2D3D_sinang;

/* layered model */

struct layer {
    double dtop; /* depth to top of layer */
    double vptop; /* P velocity at top of layer */
    double vpgrad; /* P vel gradient in layer */
    double vstop; /* S velocity at top of layer */
    double vsgrad; /* S vel gradient in layer */
    double dentop; /* density at top of layer */
    double dengrad; /* density gradient in layer */
    int plot; /* line plot style */
};

#define MAX_LAYERS 100
EXTERN_TXT struct layer model_layer[MAX_LAYERS];
EXTERN_TXT int num_layers; /* number of layers read in */
#define LAYEROFFSET 99000	/* base offset for layer model element id */


/* surface model */
/*	surface is defined by a GMT grd file */


#define SURF_REF_ABS 0
#define SURF_REF_SURF 1
#define SURF_REF_HIGH 2
#define SURF_REF_LOW 3

struct surface {
    char grd_file[MAXFILE]; /* GMT grd filename */
    struct grd_hdr *hdr; /* GMT 2-D gridfile header structure */
    float *zdata; /* z-values stored in scanline format */
    double zshift; /* amount to shift grd surface
					in km (pos = up, neg = down) */
    int iref_type; /* SURF_REF_ABS - use ref level depth
				   SURF_REF_SURF - use local depth of surface
				   SURF_REF_HIGH - use highest depth on surface
				   SURF_REF_LOW - use lowest depth on surface */
    double ref_level; /* depth from ref point to reference level
					for vel/den determination; */
    double pix_shift; /* shift in loc to account for node_offset type */
    double zmin; /* minimum z value on surface */
    double zmax; /* minimum z value on surface */
    double vptop; /* P velocity at ref level */
    double vpgrad; /* P vel gradient in layer */
    double vstop; /* S velocity at ref level */
    double vsgrad; /* S vel gradient in layer */
    double dentop; /* density at ref level */
    double dengrad; /* density gradient in layer */
    int plot; /* line plot style */
};

#define MAX_SURFACES 100
EXTERN_TXT struct surface model_surface[MAX_SURFACES];
EXTERN_TXT int num_surfaces; /* number of surfaces read in */
#define SURFACEOFFSET 3000	/* base offset for surface model element id */

/* GMT 2-D gridfile header structure */

struct grd_hdr {
    int nx; /* Number of columns */
    int ny; /* Number of rows */
    int node_offset; /* 0 for node grids, 1 for pixel grids */
    double x_min; /* Minimum x coordinate */
    double x_max; /* Maximum x coordinate */
    double y_min; /* Minimum y coordinate */
    double y_max; /* Maximum y coordinate */
    double z_min; /* Minimum z value */
    double z_max; /* Maximum z value */
    double x_inc; /* x increment */
    double y_inc; /* y increment */
    double z_scale_factor; /* grd values must be multiplied by this */
    double z_add_offset; /* After scaling, add this */
    char x_units[80]; /* units in x-direction */
    char y_units[80]; /* units in y-direction */
    char z_units[80]; /* grid value units */
    char title[80]; /* name of data set */
    char command[320]; /* name of generating command */
    char remark[160]; /* comments re this data set */

    /* Data record: float z[nx*ny] z-values stored in scanline format */

    /* Note on node_offset:
            Assume x_min = y_min = 0 and x_max = y_max = 10 and x_inc = y_inc = 1.
            For a normal node grid we have:
                    (1) nx = (x_max - x_min) / x_inc + 1 = 11
                        ny = (y_max - y_min) / y_inc + 1 = 1
                    (2) node # 0 is at (x,y) = (x_min, y_max) = (0,10) and represents the surface
                        value in a box with dimensions (1,1) centered on the node.
            For a pixel grid we have:
                    (1) nx = (x_max - x_min) / x_inc = 10
                        ny = (y_max - y_min) / y_inc = 10
                    (2) node # 0 is at (x,y) = (x_min + 0.5*x_inc, y_max - 0.5*y_inc) = (0.5, 9.5)
                        and represents the surface value in a box with dimensions (1,1)
                        centered on the node.
     */

    /*  From GMT-3_App_B.pdf:

    Scanline format means that the data are stored in rows (y = constant) going from the "top" (y = y_max (north)) to the "bottom" (y = y_min (south)). Data within each row are ordered from "left" (x = x_min (west)) to "right" (x = x_max (east)). The node_offset signals how the nodes are laid out. The grid is always defined as the intersections of all x (x = x_min, x_min + x_inc, x_min + 2*x_inc, ..., x_max) and y (y = y_min, y_min + y_inc, y_min + 2*y_inc, ..., y_max) lines.

    The two scenarios differ in which area each data point represents.

    1) Grid line registration. In this registration, the nodes are centered on the grid line intersections and the data points represent the average value in a cell of dimensions (x_inc * y_inc) centered on the nodes:

            FIGURE

    2) Pixel registration Here, the nodes are centered in the grid cells, i.e., the areas between grid lines, and the data points represent the average values within each cell:

            FIGURE

    Thus, inspecting the figures we see that in the case of grid registration the number of nodes are related to region and grid spacing by nx = (x_max - x_min) / x_inc + 1 = 4 ny = (y_max - y_min) / y_inc + 1 = 4 while for pixel registration we find nx = (x_max - x_min) / x_inc = 3 ny = (y_max - y_min) / y_inc = 3 The default registration in GMT is grid line registration. Most programs can handle both types, and for some programs like grdimage a pixel registered file makes more sense. Utility programs like grdsample and grdproject will allow you to convert from one format to the other.
     */

};


/* rough model */

#define MAX_SINUSOIDS 20
#define MAX_ROUGH 20

struct rough_bndry {
    double ref_lev; /* ref level */
    double vptop; /* P velocity at ref lev */
    double vpgrad; /* P vel gradient in layer */
    double vstop; /* S velocity at ref level */
    double vsgrad; /* S vel gradient in layer */
    double dentop; /* density at ref level */
    double dengrad; /* density gradient in layer */
    double zmax; /* max z coord of bndry */
    double zmin; /* min z coord of bndry */
    int num_sin; /* number of sin functions to define boundary */
    double amp[MAX_SINUSOIDS]; /* amplitude of sinusoid (km) */
    double wavelen[MAX_SINUSOIDS]; /* wavelength of sinusoid (km) */
    double phase[MAX_SINUSOIDS]; /* phase offset of sinusoid (km) */
    int plot; /* line plot style */
};

EXTERN_TXT struct rough_bndry model_rough[MAX_ROUGH];
EXTERN_TXT int num_rough; /* number of layers read in */

/* disk inclusions */

struct disk {
    double x; /* x loc of disk center */
    double z; /* z loc of disk center */
    double rad2; /* square of radius of disk */
    double pvel; /* P velocity in disk */
    double svel; /* S velocity in disk */
    double den; /* density in disk */
    int plot; /* line plot style */
};

#define MAX_DISKS 100
EXTERN_TXT struct disk model_disk[MAX_DISKS]; /* array of disk structures */

EXTERN_TXT int num_disks; /* number of disks read in */

/* sphere inclusions (3-D only) */

struct sphere {
    double z; /* z loc of sphere center */
    double x; /* x loc of sphere center */
    double y; /* y loc of sphere center */
    double rad2; /* square of radius of sphere */
    double pvel; /* P velocity in sphere */
    double svel; /* S velocity in sphere */
    double den; /* density in sphere */
};

/* array of sphere structures */

#define MAX_SPHERES 100
EXTERN_TXT struct sphere model_sphere[MAX_SPHERES];

EXTERN_TXT int num_spheres; /* number of spheres read in */

/** polygon model */

/* vertex structure */

struct vertex {
    struct vertex *prev; /* pointer to previous vertex */
    struct vertex *next; /* pointer to next vertex */
    int id_vtx; /* vertex identification */
    double x; /* x loc of vertex */
    double y; /* y loc of vertex */
    double z; /* z loc of vertex */
};

EXTERN_TXT struct vertex *vtx_head;
EXTERN_TXT int num_vtx; /* number of vertices read in */

/* edge structure */

struct edge {
    struct edge *prev; /* pointer to previous edge */
    struct edge *next; /* pointer to next edge */
    int id_edge; /* edge identification */
    struct vertex *v1; /* first vertex */
    struct vertex *v2; /* second vertex */
    int edge_plot; /* plotting flag, 1 = draw edge, 0 = no draw */
    int refl_mode; /* reflection mode  0 = no, 1 = yes (2D SH only) */
};

EXTERN_TXT struct edge *edge_head;
EXTERN_TXT int num_edge; /* number of edges read in */

/* 2-D or 3-D polygon structure */

struct polygon {
    struct polygon *prev; /* pointer to previous polygon */
    struct polygon *next; /* pointer to next polygon */
    int id_poly; /* polygon identification */
    int n_edges; /* number of edges */
    struct edge **edges; /* edge list */
    double xmin; /* min x val of polygon */
    double xmax; /* max x val of polygon */
    double ymin; /* min y val of polygon */
    double ymax; /* max y val of polygon */
    double zmin; /* min z val of polygon */
    double zmax; /* max z val of polygon */
    struct normals *norm; /* array of normal coefficients */
    double ref; /* velocity reference level */
    double vpref; /* P velocity at reference level */
    double vpgrad; /* P vel gradient */
    double vsref; /* S velocity at reference level */
    double vsgrad; /* S vel gradient */
    double denref; /* density at reference level */
    double dengrad; /* density gradient */
};

EXTERN_TXT struct polygon *poly_head;
EXTERN_TXT int num_poly; /* number of polygons read in */
#define POLYOFFSET 2000		/* base offset for polygon model element id */

/* 3-D solid structure */

struct solid {
    struct solid *prev; /* pointer to previous solid */
    struct solid *next; /* pointer to next solid */
    int id_solid; /* solid identification */
    int n_poly; /* number of polygonss */
    struct polygon **poly; /* polygon list */
    double xmin; /* min x val of solid */
    double xmax; /* max x val of solid */
    double ymin; /* min x val of solid */
    double ymax; /* max x val of solid */
    double zmin; /* min z val of solid */
    double zmax; /* max z val of solid */
    struct normals *norm; /* array of normal coefficients */
    double ref; /* velocity reference level */
    double vpref; /* P velocity at reference level */
    double vpgrad; /* P vel gradient */
    double vsref; /* S velocity at reference level */
    double vsgrad; /* S vel gradient */
    double denref; /* density at reference level */
    double dengrad; /* density gradient */
};

EXTERN_TXT struct solid *solid_head;
EXTERN_TXT int num_solid; /* number of solids read in */

/* edge or polygon normal coeffs (ax + by + cz + d = 0) */
struct normals {
    double a;
    double b;
    double c;
    double d;
};


/* vel fd grid model */

EXTERN_TXT int fdgrid_flag;
EXTERN_TXT int fdgrid_numx, fdgrid_numz;
EXTERN_TXT double fdgrid_xmin, fdgrid_xmax;
EXTERN_TXT double fdgrid_zmin, fdgrid_zmax;
EXTERN_TXT double fdgrid_xstep, fdgrid_zstep;
EXTERN_TXT float *fdgrid_array;
EXTERN_TXT char vfile_name[MAXFILE]; /* vgrid file name */

EXTERN_TXT double vmodel_vmean;

EXTERN_TXT double cPI; /* PI = 3.14...*/

/**************************************************************************/
/* function prototypes */
/**************************************************************************/

void set_rough_limits(struct rough_bndry *);
void read_fdiff_vel(char*);
int set_poly_normals(struct polygon *);
int set_poly_limits(struct polygon *);
int set_solid_normals(struct solid *);
int set_solid_limits(struct solid *);

int inside_poly(double, double, struct polygon *);
int inside_solid(double, double, double, struct solid *);

int read_vel_mod_input(FILE*, char*, char*, int, int);
int get_model_layer(struct layer *pm, int nlayer, char *input_line);
int get_model_disk(struct disk *pd, int ndisk, char *input_line);
double get_vel(double xpos, double ypos, double zpos, char wavetype, double *density, int idensity, int *imodel);
int get_vel_den();
double get_layer_vel(double depth, char wavetype, struct layer *pm, int nlayer, double *density, int iden, int *imodel);

double get_disk_vel(double xpos, double zpos, char wavetype, struct disk *pd, int ndisk, double *density, int iden);
double get_sphere_vel(double xpos, double ypos, double zpos, char wavetype, struct sphere *pd, int nsphere, double *density, int iden);
int get_model_surface(struct surface *, int, char *, int);
int read_grd_surface(struct surface *, int, int);
int read_grd(struct surface *ps, int imessage);
int get_model_rough(struct rough_bndry *, int, char*, FILE*);
int get_model_fdgrid(char*, FILE*);
int get_model_sphere(struct sphere *, int, char*);
int get_model_vertex(char*);
int get_model_edge(char*);
int get_model_poly(char*, FILE*);
int get_model_poly_3d(char*, FILE*);
int get_model_solid(char*, FILE*);


struct vertex *addvtx(int id_num);
//struct edge *addedge();
struct edge *addedge(int id_num);
struct polygon *addpoly(int id_num);
struct solid *addsolid(int id_num);
void disp_model_poly();
double get_poly_vel(double, double, char, double*,
        int, int *);
double get_poly_vel_2D3D(double, double, double, char,
        double*, int, int*);
double get_solid_vel(double xpos, double ypos, double zpos, char wavetype, double *density, int iden);
double get_rough_z(int nrough, double x);
double get_rough_vel(double xpos, double zpos, char wavetype, struct rough_bndry *pr, int nrough, double *density, int iden);
double get_surface_vel(double, double, double, char,
        struct surface *, int, double *, int);
double get_surface_z(int, double, double);
int free_surface(struct surface *ps);
int dump_grd(int nsurface, int idump_decimation, double x_factor, double y_factor, double z_factor, char *dump_file);

int draw_model();
int draw_dist_axis();
int draw_time_axis();

double get_fdiff_vel(double xpos, double zpos, char wavetype, double *density, int iden);



/*

EXTERN_TXT int init_graph();
EXTERN_TXT int init_disp();
EXTERN_TXT int close_graph();
EXTERN_TXT drawline();
EXTERN_TXT int draw_front();
EXTERN_TXT setline();
EXTERN_TXT setcolor();
EXTERN_TXT int set_colortable();
EXTERN_TXT gtextw();
EXTERN_TXT char keypress();
EXTERN_TXT message();
 */
