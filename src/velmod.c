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


/* 	wls Velocity Model Functions
        Functions to get velocity */

/* 	History:

        16 May 1991		original routines
        23 Sep 1991		polygon model added

 */


#define EXTERN_MODE 1

#undef EXTERN_TXT
#ifdef EXTERN_MODE
#define EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif


#include "util.h"
#include "velmod.h"


int Get2Dto3DTrans(char* input_line);

/*** function to read velocity model input */

int read_vel_mod_input(FILE* fp_input, char* param, char* line, int istat, int imessage) {

    /* read velocity model input */

    if (strcmp(param, "SURFACE") == 0) { /* read surfrace desc */
        if ((istat = get_model_surface(model_surface, num_surfaces,
                strchr(line, ' '), imessage)) < 0)
            fprintf(stderr, "ERROR: reading model surface.\n");
        else
            num_surfaces++;
    }

    if (strcmp(param, "LAYER") == 0) { /* read layer desc */
        if ((istat = get_model_layer(model_layer, num_layers,
                strchr(line, ' '))) < 0)
            fprintf(stderr, "ERROR: reading model layer.\n");
        else
            num_layers++;
    }

    if (strcmp(param, "ROUGH") == 0) { /* read rough bndry desc */
        if ((istat = get_model_rough(model_rough, num_rough,
                strchr(line, ' '), fp_input)) < 0)
            fprintf(stderr,
                "ERROR: reading model rough layer.\n");
        else
            num_rough++;
    }

    if (strcmp(param, "VGRID") == 0) { /* read vel grid desc */
        if ((istat =
                get_model_fdgrid(strchr(line, ' '), fp_input)) < 0)
            fprintf(stderr, "ERROR: reading fdgrid model.\n");
    }

    if (strcmp(param, "DISK") == 0) { /* read disk desc */
        if ((istat = get_model_disk(model_disk, num_disks,
                strchr(line, ' '))) < 0)
            fprintf(stderr, "ERROR: reading disk model.\n");
        num_disks += istat;
    }

    if (strcmp(param, "SPHERE") == 0) { /* read disk desc */
        if ((istat = get_model_sphere(model_sphere, num_spheres,
                strchr(line, ' '))) < 0)
            fprintf(stderr, "ERROR: reading sphere model.\n");
        num_spheres += istat;
    }

    if (strcmp(param, "2DTO3DTRANS") == 0) /* read 2D to 3D trans */
        if ((istat = Get2Dto3DTrans(strchr(line, ' '))) < 0)
            fprintf(stderr,
                "ERROR: reading 2D to 3D tansformation.\n");

    if (strcmp(param, "VERTEX") == 0) /* read vertex */
        if ((istat = get_model_vertex(strchr(line, ' '))) < 0)
            fprintf(stderr, "ERROR: reading vertex.\n");

    if (strcmp(param, "EDGE") == 0) /* read edge */
        if ((istat = get_model_edge(strchr(line, ' '))) < 0)
            fprintf(stderr, "ERROR: reading edge.\n");

    if (strcmp(param, "POLYGON2") == 0) /* read polygon */
        if ((istat = get_model_poly(strchr(line, ' '), fp_input)) < 0)
            fprintf(stderr, "ERROR: reading 2D polygon.\n");

    if (strcmp(param, "POLYGON3") == 0) /* read polygon */
        if ((istat = get_model_poly_3d(strchr(line, ' '), fp_input)) < 0)
            fprintf(stderr, "ERROR: reading 3D polygon.\n");

    if (strcmp(param, "SOLID") == 0) /* read polygon */
        if ((istat = get_model_solid(strchr(line, ' '), fp_input)) < 0)
            fprintf(stderr, "ERROR: reading solid.\n");

    return (istat);

}

/* function to get velocity at location */

double get_vel(double xpos, double ypos, double zpos, char wavetype, double *density, int idensity, int *imodel) {
    int iden;
    double vel;

    iden = (density != NULL && idensity) ? 1 : 0;

    *imodel = 0;

    if (num_spheres < 1 ||
            (vel = get_sphere_vel(xpos, ypos, zpos, wavetype, model_sphere,
            num_spheres, density, iden)) < 0.0) {
        if (num_disks < 1 ||
                (vel = get_disk_vel(xpos, zpos, wavetype, model_disk,
                num_disks, density, iden)) < 0.0) {
            if ((prog_mode_3d && (num_solid < 1 ||
                    (vel = get_solid_vel(xpos, ypos, zpos, wavetype,
                    density, iden)) < 0.0)) ||
                    (!prog_mode_3d && (num_poly < 1 ||
                    (vel = get_poly_vel(xpos, zpos, wavetype,
                    density, iden, imodel)) < 0.0))) {
                if (!prog_mode_Mod2D3D ||
                        (prog_mode_Mod2D3D && (num_poly < 1 ||
                        (vel = get_poly_vel_2D3D(xpos, ypos, zpos,
                        wavetype, density, iden, imodel)) < 0.0))) {
                    if (num_rough < 1 ||
                            (vel = get_rough_vel(xpos, zpos,
                            wavetype, model_rough,
                            num_rough, density, iden)) < 0.0) {
                        /*if (num_surfaces < 1 ||
                                (vel = get_surface_vel(xpos, ypos, zpos,
                                wavetype, model_surface,
                                num_surfaces, density, iden)) < 0.0)
                            {*/
                        if (fdgrid_flag < 1 ||
                                (vel = get_fdiff_vel(xpos, -zpos,
                                wavetype, density, iden)) < 0.0) {
                            vel = get_layer_vel(zpos,
                                    wavetype, model_layer,
                                    num_layers, density,
                                    iden, imodel);
                            if (vel < 0.0)
                                *imodel = -1;
                        }
                        /*}*/
                    }
                }
            }
        }
    }

    return (vel);
}

/* function to calculate velocity from gradient layer model */

double get_layer_vel(double depth, char wavetype, struct layer *pm, int nlayer, double *density, int iden, int *imodel) {
    int n;
    double vel = -1.0;

    /* reflect point above free surface */

    /*    if (depth < 0.0) */
    /*		depth = -depth; */

    /* Skip to layer containing depth */

    for (n = 0; n < nlayer - 1 && depth >= (pm + n + 1)->dtop; n++)
        ;

    /* calculate velocity */

    if (wavetype == 'P')
        vel = (pm + n)->vptop + (pm + n)->vpgrad * (depth - (pm + n)->dtop);
    else
        vel = (pm + n)->vstop + (pm + n)->vsgrad * (depth - (pm + n)->dtop);

    /* Calculate density */

    if (iden == 1)
        *density = (pm + n)->dentop + (pm + n)->dengrad * (depth - (pm + n)->dtop);

    /* set model element code */

    *imodel = LAYEROFFSET + n;

    return (vel);
}

/* function to calculate velocity from disk model */

double get_disk_vel(double xpos, double zpos, char wavetype, struct disk *pd, int ndisk, double *density, int iden) {
    int n;
    double xtemp, ztemp;
    double vel = -1.0;

    for (n = 0; n < ndisk; n++) { /* check if pos is in each disk */
        xtemp = xpos - (pd + n)->x;
        ztemp = zpos - (pd + n)->z;
        if (((xtemp * xtemp) + (ztemp * ztemp)) <= (pd + n)->rad2) {
            if (wavetype == 'P')
                vel = (pd + n)->pvel;
            else
                vel = (pd + n)->svel;
            if (iden == 1)
                *density = (pd + n)->den;
            break;
        }
    }

    return (vel);
}

/* function to calculate velocity from surface model */

double get_surface_vel(double xpos, double ypos, double zpos, char wavetype,
        struct surface *psur, int nsurface, double *density, int iden) {
    int n;
    double vel = -1.0;
    double surface_level, ref_level;
    struct surface *ps;



    /* find last surface in list above point  */

    for (n = nsurface - 1; n >= 0; n--) {
        ps = psur + n;
        if (zpos >= ps->zmin && zpos
                >= (surface_level = get_surface_z(n, xpos, ypos))) {
            if (ps->iref_type == SURF_REF_SURF)
                ref_level = surface_level + ps->ref_level;
            else
                ref_level = ps->ref_level;
            if (wavetype == 'P')
                vel = ps->vptop + ps->vpgrad * (zpos - ref_level);
            else
                vel = ps->vstop + ps->vsgrad * (zpos - ref_level);
            if (iden == 1)
                *density = ps->dentop +
                    ps->dengrad * (zpos - ref_level);
            break;
        }

    }

    /*printf("xpos %lf ypos %lf zpos %lf ns %d vel %lf\n", xpos, ypos, zpos, n, vel);*/

    return (vel);
}

/*** function to get surface boundary z level */

double get_surface_z(int nsurface, double x, double y) {

    long ioffset;
    int ixnode, iynode;
    double zdepth;
    struct surface *ps;
    struct grd_hdr *phdr;


    ps = &model_surface[nsurface];
    phdr = ps->hdr;


    /* get node coords */
    ixnode = (int) (ps->pix_shift + (x - phdr->x_min) / phdr->x_inc);
    iynode = (int) (ps->pix_shift + (y - phdr->y_min) / phdr->y_inc);
    iynode = phdr->ny - 1 - iynode; /* adjust for scanline order */

    /* check if outside grid */
    if (ixnode < 0 || ixnode >= phdr->nx || iynode < 0 || iynode >= phdr->ny)
        return (LARGE_FLOAT);

    /* get z level value at node */
    ioffset = iynode * phdr->nx + ixnode;
    zdepth = *(ps->zdata + ioffset);

    return (zdepth);
}

/* function to calculate velocity from rough layer model */

double get_rough_vel(double xpos, double zpos, char wavetype, struct rough_bndry *pr, int nrough, double *density, int iden) {
    int n;
    double vel = -1.0;


    /* check if zpos is below each bndry */

    for (n = nrough - 1; n >= 0; n--) {
        if (zpos > (pr + n)->zmin && zpos > get_rough_z(n, xpos)) {
            if (wavetype == 'P')
                vel = (pr + n)->vptop +
                    (pr + n)->vpgrad * (zpos - (pr + n)->ref_lev);
            else
                vel = (pr + n)->vstop +
                    (pr + n)->vsgrad * (zpos - (pr + n)->ref_lev);
            if (iden == 1)
                *density = (pr + n)->dentop +
                    (pr + n)->dengrad * (zpos - (pr + n)->ref_lev);
            break;
        }
    }

    return (vel);
}

/* function to calculate velocity from sphere model */

double get_sphere_vel(double xpos, double ypos, double zpos, char wavetype, struct sphere *pd, int nsphere, double *density, int iden) {
    int n;
    double xtemp, ytemp, ztemp;
    double vel = -1.0;


    for (n = 0; n < nsphere; n++) { /* check if pos is in each sphere */
        xtemp = xpos - (pd + n)->x;
        ytemp = ypos - (pd + n)->y;
        ztemp = zpos - (pd + n)->z;
        if (((xtemp * xtemp) + (ytemp * ytemp) + (ztemp * ztemp))
                <= (pd + n)->rad2) {
            if (wavetype == 'P')
                vel = (pd + n)->pvel;
            else
                vel = (pd + n)->svel;
            if (iden == 1)
                *density = (pd + n)->den;
            break;
        }
    }

    return (vel);
}

/* function to read disk inclusions  from input file */

int get_model_disk(pd, ndisk, input_line)
struct disk *pd;
int ndisk;
char *input_line;
{
    int n;
    double rad;

    n = ndisk;

    sscanf(input_line,
            "%lf %lf %lf %lf %lf %lf", &(pd + n)->x, &(pd + n)->z, &rad,
            &(pd + n)->pvel, &(pd + n)->svel, &(pd + n)->den);
    (pd + n)->rad2 = rad * rad;

    return (1);
}

/* function to read sphere inclusion  from input file */

int get_model_sphere(struct sphere *pd, int nsphere, char* input_line) {
    int n;
    double rad;

    n = nsphere;

    sscanf(input_line,
            "%lf %lf %lf %lf %lf %lf %lf", &(pd + n)->z, &(pd + n)->x, &(pd + n)->y,
            &rad, &(pd + n)->pvel, &(pd + n)->svel, &(pd + n)->den);
    (pd + n)->rad2 = rad * rad;

    printf("SPHERE  Z %lf X %lf Y %lf  rad %lf  Pvel %lf Svel %lf den %lf\n",
            (pd + n)->z, (pd + n)->x, (pd + n)->y,
            rad, (pd + n)->pvel, (pd + n)->svel, (pd + n)->den);


    return (1);
}

/* function to read 2D to 3D tansformation parameters from input file */

int Get2Dto3DTrans(char* input_line) {
    int istat;
    double angle, rad_per_deg;

    if ((istat = sscanf(input_line, "%lf %lf %lf",
            &Mod2D3D_origx, &Mod2D3D_origy, &Mod2D3D_rot)) != 3)
        return (-1);

    printf("2DTO3DTRANSFORM  X_orig %f  Y_orig %f  rot %f\n",
            Mod2D3D_origx, Mod2D3D_origy, Mod2D3D_rot);

    rad_per_deg = atan(1.0) / 45.0;
    angle = -rad_per_deg * Mod2D3D_rot;
    Mod2D3D_cosang = cos(angle);
    Mod2D3D_sinang = sin(angle);

    prog_mode_Mod2D3D = 1;



    return (1);
}


/*** function to read model surface from input file */

/*	surface is defined by a GMT grd file
        surface array is in order of reading of surfaces
        velocity at point x,y,z is determined by first surface found,
                working backwards from end of surface array, that is above point
 */

int get_model_surface(struct surface *psur, int nsurface, char *input_line, int imessage) {

    int istat;
    char str_ref_type[MAXLINE];
    struct surface *ps;


    ps = psur + nsurface;

    (ps)->plot = 1;
    if ((istat = sscanf(input_line, "%s %lf %s %lf %lf %lf %lf %lf %lf %lf %d",
            ps->grd_file, &ps->zshift, str_ref_type, &ps->ref_level,
            &ps->vptop, &ps->vpgrad,
            &ps->vstop, &ps->vsgrad,
            &ps->dentop, &ps->dengrad,
            &ps->plot))
            != 10 && istat != 11)
        return (-1);

    /* set reference level type */

    if (strcmp(str_ref_type, "REF_ABS") == 0)
        ps->iref_type = SURF_REF_ABS;
    else if (strcmp(str_ref_type, "REF_SURF") == 0)
        ps->iref_type = SURF_REF_SURF;
    else if (strcmp(str_ref_type, "REF_HIGH") == 0)
        ps->iref_type = SURF_REF_HIGH;
    else if (strcmp(str_ref_type, "REF_LOW") == 0)
        ps->iref_type = SURF_REF_LOW;
    else {
        fprintf(stderr, "ERROR: Unrecognized surface reference level type: <%s>\n",
                str_ref_type);
        return (-1);
    }


    if (read_grd_surface(ps, imessage, 1) < 0)
        return (-1);

    return (1);
}

/*** function to read GMT grd file and initialize surface */

int read_grd_surface(struct surface *ps, int imessage, int force_km) {

    int istat;
    long idatasize, npt;
    FILE *fp_grd;
    char *phline, hline[MAXLINE_LONG], filename[MAXLINE], regstr[MAXLINE];
    char *pchr, *psubstr;
    double zval;


    /* open grd file */

    if ((fp_grd = fopen(ps->grd_file, "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open surface grd file:\n");
        fprintf(stderr, "  %s\n", ps->grd_file);
        return (-1);
    }


    /* allocate and read grd file header */

    if ((ps->hdr = (struct grd_hdr *) malloc(sizeof (struct grd_hdr))) == NULL) {
        fprintf(stderr, "ERROR: Cannot allocate grd header memory.\n");
        return (-1);
    }

    if (imessage)
        printf("\nGMT grd file header:  %s\n", ps->grd_file);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    psubstr = strrchr(hline, ':');
    if (psubstr != NULL)
        strcpy((ps->hdr)->title, psubstr + 1);
    if ((pchr = strchr((ps->hdr)->title, '\n')) != NULL)
        *pchr = '\0';
    if (imessage)
        printf("\"%s\"\t\t/* Descriptive title of the data set */\n", (ps->hdr)->title);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    psubstr = strrchr(hline, ':');
    if (psubstr != NULL)
        strcpy((ps->hdr)->command, psubstr + 1);
    if ((pchr = strchr((ps->hdr)->command, '\n')) != NULL)
        *pchr = '\0';
    if (imessage)
        printf("\"%s\"\t\t/* Command line that produced the grdfile */\n",
            (ps->hdr)->command);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    psubstr = strrchr(hline, ':');
    if (psubstr != NULL)
        strcpy((ps->hdr)->remark, psubstr + 1);
    if ((pchr = strchr((ps->hdr)->remark, '\n')) != NULL)
        *pchr = '\0';
    if (imessage)
        printf("\"%s\"\t\t/* Any additional comments */\n", (ps->hdr)->remark);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    istat = sscanf(hline, "%s %s", filename, regstr);
    if (strcmp(regstr, "Normal") == 0) {
        (ps->hdr)->node_offset = 0;
        ps->pix_shift = 0.5;
    } else {
        (ps->hdr)->node_offset = 1;
        ps->pix_shift = 0.0;
    }
    if (imessage)
        printf("sscanf istat=%d\n", istat);
    printf("%d\t\t\t/* 0 for grid line reg, 1 for pixel reg */\n",
            (ps->hdr)->node_offset);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    /* skip grdfile format # line */

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    istat = sscanf(hline, "%s x_min: %lf x_max: %lf x_inc: %lf %*s %s nx: %d",
            filename, &(ps->hdr)->x_min, &(ps->hdr)->x_max, &(ps->hdr)->x_inc,
            (ps->hdr)->x_units, &(ps->hdr)->nx);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Minimum x-value of region */\n", (ps->hdr)->x_min);
        printf("%lf\t/* Maximum x-value of region */\n", (ps->hdr)->x_max);
        printf("%lf\t/* Node spacing in x-dimension */\n", (ps->hdr)->x_inc);
        printf("%s\t/* Units of the x-dimension */\n", (ps->hdr)->x_units);
        printf("%d\t\t\t/* Number of nodes in the x-dimension */\n",
                (ps->hdr)->nx);
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    istat = sscanf(hline, "%s y_min: %lf y_max: %lf y_inc: %lf %*s %s ny: %d",
            filename, &(ps->hdr)->y_min, &(ps->hdr)->y_max, &(ps->hdr)->y_inc,
            (ps->hdr)->y_units, &(ps->hdr)->ny);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Minimum y-value of region */\n", (ps->hdr)->y_min);
        printf("%lf\t/* Maximum y-value of region */\n", (ps->hdr)->y_max);
        printf("%lf\t/* Node spacing in y-dimension */\n", (ps->hdr)->y_inc);
        printf("%s\t/* Units of the y-dimension */\n", (ps->hdr)->y_units);
        printf("%d\t\t\t/* Number of nodes in the y-dimension */\n", (ps->hdr)->ny);
    }

    // check if km units  // 20161102 AJL - added
    ps->is_latlon = 1;
    if (strcmp((ps->hdr)->x_units, "km") == 0 && strcmp((ps->hdr)->y_units, "km") == 0) {
        ps->is_latlon = 0;
        printf("/* X/Y grid is kilometers. */\n");
    } else {
        printf("/* X/Y grid is lat/lon. */\n");
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    istat = sscanf(hline, "%s z_min: %lf z_max: %lf %*s %s",
            filename, &(ps->hdr)->z_min, &(ps->hdr)->z_max, (ps->hdr)->z_units);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Minimum z-value in data set */\n", (ps->hdr)->z_min);
        printf("%lf\t/* Maximum z-value in data set */\n", (ps->hdr)->z_max);
        printf("%s\t/* Units of the z-dimension */\n", (ps->hdr)->z_units);
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL)
        return (-1);
    istat = sscanf(hline, "%s  scale_factor: %lf add_offset: %lf",
            filename, &(ps->hdr)->z_scale_factor, &(ps->hdr)->z_add_offset);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Factor to multiply z-values after read */\n",
                (ps->hdr)->z_scale_factor);
        printf("%lf\t/* Offset to add to scaled z-values */\n", (ps->hdr)->z_add_offset);
    }

    /* data must be kilometers */
    if (force_km && strcmp((ps->hdr)->z_units, "km") != 0) {
        fprintf(stderr, "ERROR: Z-level data must be kilometers.\n");
        return (-1);
    }


    /* allocate and read grd data */

    idatasize = (ps->hdr)->nx * (ps->hdr)->ny;
    if ((ps->zdata = (float *) malloc((size_t) (idatasize * sizeof (float)))) == NULL) {
        fprintf(stderr, "ERROR: Cannot allocate array for grd z data.\n");
        return (-1);
    }
    for (npt = 0; npt < idatasize; npt++) {
        if (fscanf(fp_grd, " %lf", &zval) == EOF) {
            fprintf(stderr, "ERROR: Reading grd z data:\n");
            fprintf(stderr, "  %s\n", ps->grd_file);
            return (-1);
        }
        *(ps->zdata + npt) = -zval * (ps->hdr)->z_scale_factor + ps->zshift;
    }
    if (imessage)
        printf("%ld Z-level data points read.\n", idatasize);
    fclose(fp_grd);


    /* initialize surface parameters */

    if (ps->iref_type == SURF_REF_HIGH)
        ps->ref_level = -(ps->hdr)->z_max * (ps->hdr)->z_scale_factor + ps->ref_level;
    else if (ps->iref_type == SURF_REF_LOW)
        ps->ref_level = -(ps->hdr)->z_min * (ps->hdr)->z_scale_factor + ps->ref_level;
    /* set zmax and zmin in depth (pos down) */
    ps->zmax = -(ps->hdr)->z_min * (ps->hdr)->z_scale_factor + ps->zshift;
    ps->zmin = -(ps->hdr)->z_max * (ps->hdr)->z_scale_factor + ps->zshift;
    if (imessage) {
        printf("%lf\t/* Maximum depth-value of surface (after zshift) */\n", ps->zmax);
        printf("%lf\t/* Minimum depth-value of surface (after zshift)  */\n", ps->zmin);
    }

    return (0);

}

/*** function to read GMT grd file */

int read_grd(struct surface *ps, int imessage) {

    int istat;
    long idatasize;
    FILE *fp_grd;
    char *phline, hline[MAXLINE_LONG], filename[MAXLINE_LONG], regstr[MAXLINE_LONG];
    char *pchr, *psubstr;
    double zval;
    char grdfilename[MAXLINE_LONG];


    // check if binary file with ascii header file
    // grd file is GMT grd2xyz -ZTLd binary with matching ASCII GMT grdinfo header file <grd_file>.hdr
    // 20170911 AJL - added to support compact binary surface files
    //printf("DEBUG: read_grd: grid filename extension: %s\n", ps->grd_file + (strlen(ps->grd_file) - 4));
    if (strcmp(ps->grd_file + (strlen(ps->grd_file) - 4), ".bin") == 0) { // grid filename ends in ".bin", assume binary with ascii header
        ps->isHdrBinaryPair = 1;
    } else {
        ps->isHdrBinaryPair = 0;
    }


    /* open grd file */

    if (ps->isHdrBinaryPair) { // 20170911 AJL - added to support compact binary surface files
        sprintf(grdfilename, "%s.hdr", ps->grd_file);
    } else {
        sprintf(grdfilename, "%s", ps->grd_file);
    }

    if ((fp_grd = fopen(grdfilename, "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open surface grd file:\n");
        fprintf(stderr, "  %s\n", grdfilename);
        return (-1);
    }


    /* allocate and read grd file header */

    if ((ps->hdr = (struct grd_hdr *) malloc(sizeof (struct grd_hdr))) == NULL) {
        fprintf(stderr, "ERROR: Cannot allocate grd header memory.\n");
        return (-1);
    }

    if (imessage)
        printf("\nGMT grd file header:  %s\n", grdfilename);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    psubstr = strrchr(hline, ':');
    if (psubstr != NULL)
        strcpy((ps->hdr)->title, psubstr + 1);
    if ((pchr = strchr((ps->hdr)->title, '\n')) != NULL)
        *pchr = '\0';
    if (imessage)
        printf("\"%s\"\t\t/* Descriptive title of the data set */\n", (ps->hdr)->title);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    psubstr = strrchr(hline, ':');
    if (psubstr != NULL)
        strcpy((ps->hdr)->command, psubstr + 1);
    if ((pchr = strchr((ps->hdr)->command, '\n')) != NULL)
        *pchr = '\0';
    if (imessage)
        printf("\"%s\"\t\t/* Command line that produced the grdfile */\n",
            (ps->hdr)->command);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    psubstr = strrchr(hline, ':');
    if (psubstr != NULL)
        strcpy((ps->hdr)->remark, psubstr + 1);
    if ((pchr = strchr((ps->hdr)->remark, '\n')) != NULL)
        *pchr = '\0';
    if (imessage)
        printf("\"%s\"\t\t/* Any additional comments */\n", (ps->hdr)->remark);

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    istat = sscanf(hline, "%s %s", filename, regstr);
    if (strcmp(regstr, "Normal") == 0) {
        (ps->hdr)->node_offset = 0;
        ps->pix_shift = 0.5;
    } else {
        (ps->hdr)->node_offset = 1;
        ps->pix_shift = 0.0;
    }
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%d\t\t\t/* 0 for grid line reg, 1 for pixel reg */\n",
                (ps->hdr)->node_offset);
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    /* skip grdfile format # line */

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    istat = sscanf(hline, "%s x_min: %lf x_max: %lf x_inc: %lf %*s %s nx: %d",
            filename, &(ps->hdr)->x_min, &(ps->hdr)->x_max, &(ps->hdr)->x_inc,
            (ps->hdr)->x_units, &(ps->hdr)->nx);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Minimum x-value of region */\n", (ps->hdr)->x_min);
        printf("%lf\t/* Maximum x-value of region */\n", (ps->hdr)->x_max);
        printf("%lf\t/* Node spacing in x-dimension */\n", (ps->hdr)->x_inc);
        printf("%s\t/* Units of the x-dimension */\n", (ps->hdr)->x_units);
        printf("%d\t\t\t/* Number of nodes in the x-dimension */\n", (ps->hdr)->nx);
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    istat = sscanf(hline, "%s y_min: %lf y_max: %lf y_inc: %lf %*s %s ny: %d",
            filename, &(ps->hdr)->y_min, &(ps->hdr)->y_max, &(ps->hdr)->y_inc,
            (ps->hdr)->y_units, &(ps->hdr)->ny);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Minimum y-value of region */\n", (ps->hdr)->y_min);
        printf("%lf\t/* Maximum y-value of region */\n", (ps->hdr)->y_max);
        printf("%lf\t/* Node spacing in y-dimension */\n", (ps->hdr)->y_inc);
        printf("%s\t/* Units of the y-dimension */\n", (ps->hdr)->y_units);
        printf("%d\t\t\t/* Number of nodes in the y-dimension */\n", (ps->hdr)->ny);
    }

    // check if km units  // 20161102 AJL - added
    ps->is_latlon = 1;
    if (strcmp((ps->hdr)->x_units, "km") == 0 && strcmp((ps->hdr)->y_units, "km") == 0) {
        ps->is_latlon = 0;
        if (imessage)
            printf("/* X/Y grid is kilometers. */\n");
    } else {
        if (imessage)
            printf("/* X/Y grid is lat/lon. */\n");
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    istat = sscanf(hline, "%s z_min: %lf z_max: %lf %*s %s",
            filename, &(ps->hdr)->z_min, &(ps->hdr)->z_max, (ps->hdr)->z_units);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Minimum z-value in data set */\n", (ps->hdr)->z_min);
        printf("%lf\t/* Maximum z-value in data set */\n", (ps->hdr)->z_max);
        printf("%s\t/* Units of the z-dimension */\n", (ps->hdr)->z_units);
    }

    if ((phline = fgets(hline, MAXLINE_LONG, fp_grd)) == NULL) {
        fclose(fp_grd);
        return (-1);
    }
    if (imessage)
        printf("phline: \"%s\"\n", phline);
    istat = sscanf(hline, "%s  scale_factor: %lf add_offset: %lf",
            filename, &(ps->hdr)->z_scale_factor, &(ps->hdr)->z_add_offset);
    if (imessage) {
        printf("sscanf istat=%d\n", istat);
        printf("%lf\t/* Factor to multiply z-values after read */\n",
                (ps->hdr)->z_scale_factor);
        printf("%lf\t/* Offset to add to scaled z-values */\n", (ps->hdr)->z_add_offset);
    }


    // allocate and read grd data

    idatasize = (ps->hdr)->nx * (ps->hdr)->ny;
    if ((ps->zdata = (float *) malloc((size_t) (idatasize * sizeof (float)))) == NULL) {
        fprintf(stderr, "ERROR: Cannot allocate array for grd z data.\n");
        fclose(fp_grd);
        return (-1);
    }

    if (ps->isHdrBinaryPair) { // 20170911 AJL - added to support compact binary surface files
        fclose(fp_grd);
        sprintf(grdfilename, "%s", ps->grd_file);
        if ((fp_grd = fopen(grdfilename, "r")) == NULL) {
            fprintf(stderr, "ERROR: Cannot open binary surface grd file:\n");
            fprintf(stderr, "  %s\n", grdfilename);
            return (-1);
        }
        long npt_read = 0;
        while (npt_read < idatasize) {
            if (fread(&zval, sizeof (double), 1, fp_grd) != 1) {
                fprintf(stderr, "ERROR: Reading binary grd z data at npt %ld/%ld:\n", npt_read, idatasize);
                fprintf(stderr, "  %s\n", grdfilename);
                fclose(fp_grd);
                return (-1);
            }
            //if (npt_read < 50) printf("DEBUG: Z %f -> %f\n", zval, zval * (ps->hdr)->z_scale_factor + ps->zshift);
            *(ps->zdata + npt_read) = zval * (ps->hdr)->z_scale_factor + ps->zshift;
            npt_read++;
        }
        if (imessage)
            printf("%ld/%ld Z-level binary data points read.\n", npt_read, idatasize);
    } else { // single GMT ASCII file
        long npt_read = 0;
        char cchr;
        while (npt_read < idatasize) {
            if ((istat = fscanf(fp_grd, "%lf ", &zval)) == 0) { // 20170807 AJL - Added for case of additional header lines
                //printf("istat=%d\n", istat);
                // assume more header lines
                //fgets(hline, MAXLINE_LONG, fp_grd);
                while ((cchr = fgetc(fp_grd)) != '\n' && cchr != EOF) {
                    ;
                }
                continue;
            }
            if (istat == EOF) {
                fprintf(stderr, "ERROR: Reading ASCII grd z data at npt %ld/%ld:\n", npt_read, idatasize);
                fprintf(stderr, "  %s\n", grdfilename);
                fclose(fp_grd);
                return (-1);
            }
            //if (npt_read < 50) printf("DEBUG: Z %f -> %f\n", zval, zval * (ps->hdr)->z_scale_factor + ps->zshift);
            *(ps->zdata + npt_read) = zval * (ps->hdr)->z_scale_factor + ps->zshift;
            npt_read++;
        }
        if (imessage)
            printf("%ld/%ld Z-level ASCII data points read.\n", npt_read, idatasize);
    }
    fclose(fp_grd);


    return (0);

}

/*** function to free surface structure */

int free_surface(struct surface *ps) {

    if (ps->hdr != NULL)
        free(ps->hdr);
    if (ps->zdata != NULL)
        free(ps->zdata);

    return (0);
}

/*** function to write surace data to xyz file */

int dump_grd(int nsurface, int idump_decimation, double x_factor, double y_factor, double z_factor, char *dump_file) {

    double x, y;
    int npoints;
    float fdata[4], value, fmax;
    FILE *fp_out;
    struct surface *ps;


    ps = &model_surface[nsurface];


    if ((fp_out = fopen(dump_file, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open dump file:\n");
        fprintf(stderr, "  %s\n", dump_file);
        return (-1);
    }

    // skip header record (used later to store number of samples taken)
    fseek(fp_out, 4 * sizeof (float), SEEK_SET);

    // write data
    fmax = -LARGE_FLOAT;
    npoints = 0;
    for (x = (ps->hdr)->x_min; x < (ps->hdr)->x_max; x += (double) idump_decimation * (ps->hdr)->x_inc) {
        fdata[0] = x_factor * x;
        for (y = (ps->hdr)->y_min; y < (ps->hdr)->y_max; y += (double) idump_decimation * (ps->hdr)->y_inc) {
            fdata[1] = y_factor * y;
            value = z_factor * get_surface_z(nsurface, x, y);
            fdata[2] = fdata[3] = value;
            fwrite(fdata, sizeof (float), 4, fp_out);
            if (value > fmax)
                fmax = value;
            npoints++;
        }
    }

    // write header
    fseek(fp_out, 0, SEEK_SET);
    fwrite(&npoints, sizeof (int), 1, fp_out);
    fwrite(&fmax, sizeof (float), 1, fp_out);

    fclose(fp_out);

    return (0);

}

/*** function to read GMT grd file */

int convGridTokm(struct surface *ps, int imessage) {
    double vtmp;

    if (strcmp((ps->hdr)->x_units, "km") == 0) {
        ;
    } else if (strcmp((ps->hdr)->x_units, "meters") == 0) {
        (ps->hdr)->x_min /= 1000.0;
        (ps->hdr)->x_max /= 1000.0;
        (ps->hdr)->x_inc /= 1000.0;
        strcpy((ps->hdr)->x_units, "km");
    } else {
        fprintf(stderr, "ERROR: unrecognized grid x units: %s.\n",
                (ps->hdr)->x_units);
        return (-1);
    }

    if (strcmp((ps->hdr)->y_units, "km") == 0) {
        ;
    } else if (strcmp((ps->hdr)->y_units, "meters") == 0) {
        (ps->hdr)->y_min /= 1000.0;
        (ps->hdr)->y_max /= 1000.0;
        (ps->hdr)->y_inc /= 1000.0;
        strcpy((ps->hdr)->y_units, "km");
    } else {
        fprintf(stderr, "ERROR: unrecognized grid y units: %s.\n",
                (ps->hdr)->y_units);
        return (-1);
    }

    //invert y limits
    vtmp = (ps->hdr)->y_min;
    (ps->hdr)->y_min = -(ps->hdr)->y_max;
    (ps->hdr)->y_max = -vtmp;

    return (0);

}

/*** function to read model layer from input file */

int get_model_layer(struct layer *pm, int nlayer, char *input_line) {
    int n;

    n = nlayer;

    (pm + n)->plot = 1;
    sscanf(input_line, "%lf %lf %lf %lf %lf %lf %lf %d",
            &(pm + n)->dtop, &(pm + n)->vptop, &(pm + n)->vpgrad,
            &(pm + n)->vstop, &(pm + n)->vsgrad, &(pm + n)->dentop, &(pm + n)->dengrad,
            &(pm + n)->plot);

    return (1);
}

/*** function to read model rough layer from input file */

int get_model_rough(struct rough_bndry *pm, int nrough,
        char* input_line, FILE* fp)
/*!struct rough_bndry *pm;
int nrough;
char *input_line;
FILE *fp;*/ {
    int n, nsin;

    n = nrough;

    /* read layer description */

    (pm + n)->plot = 1;
    sscanf(input_line, "%lf %lf %lf %lf %lf %lf %lf %d %d",
            &(pm + n)->ref_lev, &(pm + n)->vptop, &(pm + n)->vpgrad,
            &(pm + n)->vstop, &(pm + n)->vsgrad, &(pm + n)->dentop, &(pm + n)->dengrad,
            &(pm + n)->num_sin, &(pm + n)->plot);

    /* read each sinusoid description */

    for (nsin = 0; nsin < (pm + n)->num_sin; nsin++)
        if ((fscanf(fp, "%lf %lf %lf", &(pm + n)->amp[nsin],
                &(pm + n)->wavelen[nsin], &(pm + n)->phase[nsin])) != 3)
            return (-1);

    return (1);
}

/*** function to set approx max and min range of rough boundary */

void set_rough_limits(struct rough_bndry *pr) {
    int nrough;
    double x, xmin, xmax, dx;
    double roughmin, roughmax, ztest;

    for (nrough = 0; nrough < num_rough; nrough++) {

        xmin = min_x_cut;
        xmax = max_x_cut;
        dx = (xmax - xmin) / 10000.0;
        roughmin = roughmax = (pr + nrough)->ref_lev;

        for (x = xmin; x < xmax + dx; x += dx) {
            if ((ztest = get_rough_z(nrough, x)) < roughmin)
                roughmin = ztest;
            if (ztest > roughmax)
                roughmax = ztest;
        }

        (pr + nrough)->zmin = roughmin;
        (pr + nrough)->zmax = roughmax;

    }
}

/*** function to get rough boundary z level */

double get_rough_z(nrough, x)
int nrough;
double x;
{
    int nsin;
    double z;
    struct rough_bndry *rp;

    rp = &model_rough[nrough];

    z = rp->ref_lev;

    for (nsin = 0; nsin < rp->num_sin; nsin++) {
        z += (rp->amp[nsin] / 2.0) *
                sin(2.0 * cPI * (x - rp->phase[nsin]) / rp->wavelen[nsin]);
    }

    return (z);
}

/** function to read vel grid file params fom input line */

int get_model_fdgrid(char* in_line, FILE* fp) {

    /* read description */

    sscanf(in_line, "%d %lf %lf %lf %lf %s",
            &fdgrid_flag,
            &fdgrid_xmin, &fdgrid_xmax, &fdgrid_zmin, &fdgrid_zmax,
            vfile_name);

    min_x_cut = fdgrid_xmin;
    max_x_cut = fdgrid_xmax;
    min_y_cut = -1.0;
    max_y_cut = 1.0;
    min_z_cut = fdgrid_zmin;
    max_z_cut = fdgrid_zmax;


    /* read fdgrid file */

    read_fdiff_vel(vfile_name);

    return (1);
}

/*** function to read model vertices from input file */

int get_model_vertex(char* input_line) {
    int ident;
    double x, y, z;
    struct vertex *vtxaddr;

    if ((sscanf(input_line, "%d %lf %lf %lf", &ident, &z, &x, &y)) != 4)
        return (-1);

    if ((vtxaddr = addvtx(ident)) == NULL)
        return (-2);
    num_vtx++;
    vtxaddr->id_vtx = ident;
    vtxaddr->x = x;
    vtxaddr->y = y;
    /*	vtxaddr->z = -z; */
    vtxaddr->z = z;

    return (1);
}

/*** function to read model edge from input file */

int get_model_edge(char* input_line) {
    int ident, numarg, iplot, ireflmode;
    int vtx1, vtx2;
    struct edge *edgeaddr;
    struct vertex *vtxaddr;

    if ((numarg =
            sscanf(input_line, "%d %d %d %d %d",
            &ident, &vtx1, &vtx2, &iplot, &ireflmode)) < 3)
        return (-1);

    if (numarg < 4)
        iplot = 1;

    if (numarg < 5)
        ireflmode = 0;

    if ((edgeaddr = addedge(ident)) == NULL)
        return (-2);
    num_edge++;
    edgeaddr->id_edge = ident;
    edgeaddr->edge_plot = iplot;
    edgeaddr->refl_mode = ireflmode;

    /* find addresses of edge vertices */

    edgeaddr->v1 = edgeaddr->v2 = NULL;
    vtxaddr = vtx_head;
    do {
        if (vtxaddr->id_vtx == vtx1)
            edgeaddr->v1 = vtxaddr;
        if (vtxaddr->id_vtx == vtx2)
            edgeaddr->v2 = vtxaddr;
        vtxaddr = vtxaddr->next;
    } while (vtxaddr != vtx_head &&
            (edgeaddr->v1 == NULL || edgeaddr->v2 == NULL));

    if (edgeaddr->v1 == NULL)
        fprintf(stderr, "ERROR: cannot find vertex %d.\n", vtx1);
    if (edgeaddr->v2 == NULL)
        fprintf(stderr, "ERROR: cannot find vertex %d.\n", vtx2);

    return (1);
}

/*** function to read model 2-D polygon from input file */

int get_model_poly(char* input_line, FILE* fp) {
    int ident, num_edges, nedge, id_edge;
    double ref, vpref, vpgrad, vsref, vsgrad, denref, dengrad;
    struct edge *edgeaddr;
    struct polygon *polyaddr;

    if ((sscanf(input_line, "%d %d %lf %lf %lf %lf %lf %lf %lf",
            &ident, &num_edges,
            &ref, &vpref, &vpgrad, &vsref, &vsgrad, &denref, &dengrad)) != 9)
        return (-1);

    /* skip rest of input and return if in 3D mode */

    if (prog_mode_3d && !prog_mode_Mod2D3D) {
        for (nedge = 0;
                nedge < num_edges && fscanf(fp, "%d", &id_edge); nedge++);
        return (0);
    }

    if (num_edges < 3)
        fprintf(stderr,
            "Warning polygon %d has less than 3 edges!\n", ident);


    /* add polygon structure to polygon list */

    if ((polyaddr = addpoly(ident)) == NULL)
        return (-2);
    num_poly++;
    polyaddr->id_poly = ident;
    polyaddr->n_edges = num_edges;
    polyaddr->ref = ref;
    polyaddr->vpref = vpref;
    polyaddr->vpgrad = vpgrad;
    polyaddr->vsref = vsref;
    polyaddr->vsgrad = vsgrad;
    polyaddr->denref = denref;
    polyaddr->dengrad = dengrad;

    /* allocate space for edge pointer array */

    if ((polyaddr->edges = (struct edge **)
            malloc((size_t) num_edges * sizeof (struct edge *))) == NULL)
        return (-3);

    /* read and find addresses of edges */

    for (nedge = 0; nedge < num_edges; nedge++) {
        if ((fscanf(fp, "%d", &id_edge)) != 1)
            return (-4);
        polyaddr->edges[nedge] = NULL;
        edgeaddr = edge_head;
        do {
            if (edgeaddr->id_edge == id_edge)
                polyaddr->edges[nedge] = edgeaddr;
            edgeaddr = edgeaddr->next;
        } while (edgeaddr != edge_head && polyaddr->edges[nedge] == NULL);

        if (polyaddr->edges[nedge] == NULL)
            fprintf(stderr, "ERROR: cannot find edge %d.\n", id_edge);
    }

    /* set edge normals */

    set_poly_normals(polyaddr);

    /* set polygon x, y and z limits */

    set_poly_limits(polyaddr);

    return (1);
}

/*** function to read model 3-D polygon from input file */

int get_model_poly_3d(char* input_line, FILE* fp) {
    int ident, num_edges, nedge, id_edge;
    struct edge *edgeaddr;
    struct polygon *polyaddr;

    if ((sscanf(input_line, "%d %d",
            &ident, &num_edges)) != 2)
        return (-1);

    /* skip rest of input and return if in 2D mode */

    if (!prog_mode_3d) {
        for (nedge = 0;
                nedge < num_edges && fscanf(fp, "%d", &id_edge); nedge++);
        return (0);
    }

    if (num_edges < 3)
        fprintf(stderr,
            "Warning polygon %d has less than 3 edges!\n", ident);

    /* add polygon structure to polygon list */

    if ((polyaddr = addpoly(ident)) == NULL)
        return (-2);
    num_poly++;
    polyaddr->id_poly = ident;
    polyaddr->n_edges = num_edges;

    /* allocate space for edge pointer array */

    if ((polyaddr->edges = (struct edge **)
            malloc((size_t) num_edges * sizeof (struct edge *))) == NULL)
        return (-3);

    /* read and find addresses of edges */

    for (nedge = 0; nedge < num_edges; nedge++) {
        if ((fscanf(fp, "%d", &id_edge)) != 1)
            return (-4);
        polyaddr->edges[nedge] = NULL;
        edgeaddr = edge_head;
        do {
            if (edgeaddr->id_edge == id_edge)
                polyaddr->edges[nedge] = edgeaddr;
            edgeaddr = edgeaddr->next;
        } while (edgeaddr != edge_head && polyaddr->edges[nedge] == NULL);

        if (polyaddr->edges[nedge] == NULL)
            fprintf(stderr, "ERROR: cannot find edge %d.\n", id_edge);
    }

    /* set polygon x, y and z limits */

    set_poly_limits(polyaddr);

    return (1);
}

/*** function to read model 3-D solid from input file */

int get_model_solid(char* input_line, FILE* fp) {
    int ident, num_poly, npoly, id_poly;
    double ref, vpref, vpgrad, vsref, vsgrad, denref, dengrad;
    struct polygon *polyaddr;
    struct solid *solidaddr;

    if ((sscanf(input_line, "%d %d %lf %lf %lf %lf %lf %lf %lf",
            &ident, &num_poly,
            &ref, &vpref, &vpgrad, &vsref, &vsgrad, &denref, &dengrad)) != 9)
        return (-1);

    /* skip rest of input and return if in 2D mode */

    if (!prog_mode_3d) {
        for (npoly = 0;
                npoly < num_poly && fscanf(fp, "%d", &id_poly); npoly++);
        return (0);
    }

    if (num_poly < 4)
        fprintf(stderr,
            "Warning solid %d has less than 3 polygons!\n", ident);


    /* add solid structure to solid list */

    if ((solidaddr = addsolid(ident)) == NULL)
        return (-2);
    num_solid++;
    solidaddr->id_solid = ident;
    solidaddr->n_poly = num_poly;
    solidaddr->ref = ref;
    solidaddr->vpref = vpref;
    solidaddr->vpgrad = vpgrad;
    solidaddr->vsref = vsref;
    solidaddr->vsgrad = vsgrad;
    solidaddr->denref = denref;
    solidaddr->dengrad = dengrad;

    /* allocate space for polygon pointer array */

    if ((solidaddr->poly = (struct polygon **)
            malloc((size_t) num_poly * sizeof (struct polygon *))) == NULL)
        return (-3);

    /* read and find addresses of polygons */

    for (npoly = 0; npoly < num_poly; npoly++) {
        if ((fscanf(fp, "%d", &id_poly)) != 1)
            return (-4);
        solidaddr->poly[npoly] = NULL;
        polyaddr = poly_head;
        do {
            if (polyaddr->id_poly == id_poly)
                solidaddr->poly[npoly] = polyaddr;
            polyaddr = polyaddr->next;
        } while (polyaddr != poly_head && solidaddr->poly[npoly] == NULL);

        if (solidaddr->poly[npoly] == NULL)
            fprintf(stderr, "ERROR: cannot find poly %d.\n", id_poly);
    }

    /* set polygon normals */

    set_solid_normals(solidaddr);

    /* set solid x, y and z limits */

    set_solid_limits(solidaddr);

    return (1);
}


/*** function to add a new element to vertex list */

/* (double-thread, circular linked list) */

struct vertex *
addvtx(id_num)
int id_num;
{

    struct vertex *addr;

    addr = (struct vertex *) malloc(sizeof (struct vertex));
    if (addr == NULL) {
        nll_puterr("ERROR: adding vertex, no memory?");
        return (addr);
    }

    if (vtx_head == NULL) { /* no elements in list */
        vtx_head = addr->next = addr->prev = addr;
    } else { /* add new element to end of list */
        addr->prev = vtx_head->prev;
        addr->next = vtx_head;
        addr->next->prev = addr;
        addr->prev->next = addr;
    }

    return (addr);
}


/*** function to add a new element to edge list */

/* (double-thread, circular linked list) */

struct edge *
addedge(id_num)
int id_num;
{

    struct edge *addr;

    addr = (struct edge *) malloc(sizeof (struct edge));
    if (addr == NULL) {
        nll_puterr("ERROR: adding edge, no memory?");
        return (addr);
    }

    if (edge_head == NULL) { /* no elements in list */
        edge_head = addr->next = addr->prev = addr;
    } else { /* add new element to end of list */
        addr->prev = edge_head->prev;
        addr->next = edge_head;
        addr->next->prev = addr;
        addr->prev->next = addr;
    }

    return (addr);
}


/*** function to add a new element to polygon list */

/* (double-thread, circular linked list) */

struct polygon *
addpoly(id_num)
int id_num;
{

    struct polygon *addr;

    addr = (struct polygon *) malloc(sizeof (struct polygon));
    if (addr == NULL) {
        nll_puterr("ERROR: adding polygon, no memory?");
        return (addr);
    }

    if (poly_head == NULL) { /* no elements in list */
        poly_head = addr->next = addr->prev = addr;
    } else { /* add new element to end of list */
        addr->prev = poly_head->prev;
        addr->next = poly_head;
        addr->next->prev = addr;
        addr->prev->next = addr;
    }

    return (addr);
}


/*** function to add a new element to solid list */

/* (double-thread, circular linked list) */

struct solid *
addsolid(id_num)
int id_num;
{

    struct solid *addr;

    addr = (struct solid *) malloc(sizeof (struct solid));
    if (addr == NULL) {
        nll_puterr("ERROR: adding solid, no memory?");
        return (addr);
    }

    if (solid_head == NULL) { /* no elements in list */
        solid_head = addr->next = addr->prev = addr;
    } else { /* add new element to end of list */
        addr->prev = solid_head->prev;
        addr->next = solid_head;
        addr->next->prev = addr;
        addr->prev->next = addr;
    }

    return (addr);
}

/*** function to display model vertices */

int disp_vertices() {
    struct vertex *addr;

    if ((addr = vtx_head) == NULL) {
        printf("VERTEX  No vertices read.\n");
        return (0);
    }

    printf("VERTICES (%d read)\n  ", num_vtx);
    do {
        printf("v%d:{%.2f,%.2f}  ", addr->id_vtx, addr->x, addr->z);
    } while ((addr = addr->next) != vtx_head);
    printf("\n");

    return (0);

}

/*** function to display model edges */

int disp_edges() {
    struct edge *addr;

    if ((addr = edge_head) == NULL) {
        printf("EDGE  No edges read.\n");
        return (0);
    }

    printf("EDGES (%d read)\n  ", num_edge);
    do {
        printf("e%d:{v%d,v%d}p%d,r%d  ",
                addr->id_edge, addr->v1->id_vtx, addr->v2->id_vtx,
                addr->edge_plot, addr->refl_mode);
    } while ((addr = addr->next) != edge_head);
    printf("\n");

    return (0);

}

/*** function to display model polygons */

int disp_polygons() {
    int nedge;
    struct polygon *addr;

    if ((addr = poly_head) == NULL) {
        printf("   POLYGON  No polygons read.\n");
        return (0);
    }

    printf("POLYGON (%d read)\n", num_poly);

    do {
        printf("  %d : ", addr->id_poly);
        if (!prog_mode_3d) {
            printf("ref_level=%5.2lf\n", addr->ref);
            printf("    Vp=%5.3lf dV=%6.4lf Vs=%5.3lf dV=%6.4lf",
                    addr->vpref, addr->vpgrad, addr->vsref, addr->vsgrad);
            printf(" Den=%5.2lf dDen=%6.4lf\n", addr->denref, addr->dengrad);
        }
        printf("    {");
        for (nedge = 0; nedge < addr->n_edges; nedge++) {
            if (nedge > 0)
                printf(", ");
            printf("e%d", addr->edges[nedge]->id_edge);
        }
        printf("}\n");
    } while ((addr = addr->next) != poly_head);

    return (0);

}

/*** function to display model solids */

int disp_solids() {
    int npoly;
    struct solid *addr;

    if ((addr = solid_head) == NULL) {
        printf("   SOLID  No solids read.\n");
        return (0);
    }

    printf("SOLID (%d read)\n", num_solid);

    do {
        printf("  %d : ", addr->id_solid);
        printf("ref_level=%5.2lf\n", addr->ref);
        printf("    Vp=%5.3lf dV=%6.4lf Vs=%5.3lf dV=%6.4lf",
                addr->vpref, addr->vpgrad, addr->vsref, addr->vsgrad);
        printf(" Den=%5.2lf dDen=%6.4lf\n", addr->denref, addr->dengrad);
        printf("    {");
        for (npoly = 0; npoly < addr->n_poly; npoly++) {
            if (npoly > 0)
                printf(", ");
            printf("p%d", addr->poly[npoly]->id_poly);
        }
        printf("}\n");
    } while ((addr = addr->next) != solid_head);

    return (0);

}


/*** function to set normals to edges */
/*	  Normals are used to determine if a point is inside polygon */

/*		if eq. of line is F(x,y) = Ax - By + C = 0, then point is to one side
                of the line or the other if F(x,y) > or < 0.
                This function determines A, B and C so that
                points inside the polygon will satisfy F(x,y) >= 0.
 */

int set_poly_normals(struct polygon *addr) {
    int nedge, nedge_test;
    double x1, z1, x2, z2, xtest, ztest;
    double a, b, c;
    struct vertex *vaddr1, *vaddr2, *vaddr_test;

    /* allocate space for normal coeff structures */

    if ((addr->norm = (struct normals *)
            malloc(addr->n_edges * sizeof (struct normals))) == NULL)
        nll_puterr("ERROR: allocating normals memory.");

    /* set normals for each edge */

    for (nedge = 0; nedge < addr->n_edges; nedge++) {

        /* determine coefficents */

        vaddr1 = addr->edges[nedge]->v1;
        x1 = vaddr1->x;
        z1 = vaddr1->z;
        vaddr2 = addr->edges[nedge]->v2;
        x2 = vaddr2->x;
        z2 = vaddr2->z;
        if (x1 == x2) {
            a = 1.0;
            b = 0.0;
            c = -x1;
        } else {
            a = (z1 - z2) / (x1 - x2);
            b = -1;
            c = z1 - a * x1;
        }

        /* find test vertex not on current edge */

        if ((nedge_test = nedge + addr->n_edges / 2) >= addr->n_edges)
            nedge_test = nedge_test % addr->n_edges;
        vaddr_test = addr->edges[nedge_test]->v1;
        if ((vaddr1->id_vtx == vaddr_test->id_vtx) ||
                (vaddr2->id_vtx == vaddr_test->id_vtx))
            vaddr_test = addr->edges[nedge_test]->v2;
        xtest = vaddr_test->x;
        ztest = vaddr_test->z;

        /* adjust coeff so test vertex gives > 0 in line eq. */

        if ((a * xtest + b * ztest + c) < 0.0) {
            a = -a;
            b = -b;
            c = -c;
        }

        /* store coeff */

        addr->norm[nedge].a = a;
        addr->norm[nedge].b = b;
        addr->norm[nedge].c = c;

        printf("EDGE NORM p%d e%d  %f %f %f\n", addr->id_poly,
                addr->edges[nedge]->id_edge, a, b, c);

    }

    return (0);


}

/*** function to set normals to polygons */
/*	  Normals are used to determine if a point is inside solid */

/*		if eq. of polygon is F(x,y) = Ax + By + Cz + D = 0, then point
                is to one side of the line or the other if F(x,y) > or < 0.
                This function determines A, B, C and D so that
                points inside the solid will satisfy F(x,y) >= 0.
 */

int set_solid_normals(struct solid *addr) {
    int npoly, npoly_test;
    int nedge_test, nedge_curr;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, xtest, ytest, ztest;
    double a, b, c, d;
    struct vertex *vaddr1, *vaddr2, *vaddr3, *vaddr_test = NULL, *vaddr_curr;

    /* allocate space for normal coeff structures */

    if ((addr->norm = (struct normals *)
            malloc(addr->n_poly * sizeof (struct normals))) == NULL)
        nll_puterr("ERROR: allocating normals memory.");

    /* set normals for each polygon */

    for (npoly = 0; npoly < addr->n_poly; npoly++) {

        /* select 3 vertexes from polygopn edges */

        vaddr1 = addr->poly[npoly]->edges[1]->v1;
        vaddr2 = addr->poly[npoly]->edges[1]->v2;
        vaddr3 = addr->poly[npoly]->edges[2]->v1;
        if (vaddr3 == vaddr2 || vaddr3 == vaddr1)
            vaddr3 = addr->poly[npoly]->edges[2]->v2;

        /* determine coefficents */

        x1 = vaddr1->x;
        y1 = vaddr1->y;
        z1 = vaddr1->z;
        x2 = vaddr2->x;
        y2 = vaddr2->y;
        z2 = vaddr2->z;
        x3 = vaddr3->x;
        y3 = vaddr3->y;
        z3 = vaddr3->z;

        a = y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2);
        b = z1 * (x2 - x3) + z2 * (x3 - x1) + z3 * (x1 - x2);
        c = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
        d = -x1 * (y2 * z3 - y3 * z2)
                - x2 * (y3 * z1 - y1 * z3)
                - x3 * (y1 * z2 - y2 * z1);

        /* find test vertex not on current polygon */

        for (npoly_test = 0; npoly_test < addr->n_poly; npoly_test++) {

            for (nedge_test = 0;
                    nedge_test < addr->poly[npoly_test]->n_edges; nedge_test++) {

                vaddr_test = addr->poly[npoly_test]->edges[nedge_test]->v1;

                for (nedge_curr = 0;
                        nedge_curr < addr->poly[npoly]->n_edges; nedge_curr++) {

                    vaddr_curr = addr->poly[npoly]->edges[nedge_curr]->v1;
                    if (vaddr_curr->id_vtx == vaddr_test->id_vtx)
                        goto fail1;
                    vaddr_curr = addr->poly[npoly]->edges[nedge_curr]->v2;
                    if (vaddr_curr->id_vtx == vaddr_test->id_vtx)
                        goto fail1;
                }
                goto test_exit;
fail1:

                vaddr_test = addr->poly[npoly_test]->edges[nedge_test]->v2;

                for (nedge_curr = 0;
                        nedge_curr < addr->poly[npoly]->n_edges; nedge_curr++) {

                    vaddr_curr = addr->poly[npoly]->edges[nedge_curr]->v1;
                    if (vaddr_curr->id_vtx == vaddr_test->id_vtx)
                        goto fail2;
                    vaddr_curr = addr->poly[npoly]->edges[nedge_curr]->v2;
                    if (vaddr_curr->id_vtx == vaddr_test->id_vtx)
                        goto fail2;
                }
                goto test_exit;
fail2:
                continue;
            }
        }

        fprintf(stderr, "ERROR: test vertex not found.");

test_exit:

        xtest = vaddr_test->x;
        ytest = vaddr_test->y;
        ztest = vaddr_test->z;

        /* adjust coeff so test vertex gives > 0 in line eq. */

        if ((a * xtest + b * ytest + c * ztest + d) < 0.0) {
            a = -a;
            b = -b;
            c = -c;
            d = -d;
        }

        /* store coeff */

        addr->norm[npoly].a = a;
        addr->norm[npoly].b = b;
        addr->norm[npoly].c = c;
        addr->norm[npoly].d = d;

        /*		printf("POLY NORM s%d p%d  %f %f %f %f\n", addr->id_solid,
                                addr->poly[npoly]->id_poly, a, b, c, d); */

    }

    return (0);


}

/*** function to set polygon limits */

int set_poly_limits(struct polygon *addr) {
    int nedge;
    double xmin, ymin, zmin, xmax, ymax, zmax, temp;

    xmin = ymin = zmin = REAL_MAX;
    xmax = ymax = zmax = -REAL_MAX;

    /* check limits for each vertex of each edge */

    for (nedge = 0; nedge < addr->n_edges; nedge++) {
        if ((temp = addr->edges[nedge]->v1->x) > xmax)
            xmax = temp;
        if (temp < xmin)
            xmin = temp;
        if ((temp = addr->edges[nedge]->v2->x) > xmax)
            xmax = temp;
        if (temp < xmin)
            xmin = temp;
        if ((temp = addr->edges[nedge]->v1->y) > ymax)
            ymax = temp;
        if (temp < ymin)
            ymin = temp;
        if ((temp = addr->edges[nedge]->v2->y) > ymax)
            ymax = temp;
        if (temp < ymin)
            ymin = temp;
        if ((temp = addr->edges[nedge]->v1->z) > zmax)
            zmax = temp;
        if (temp < zmin)
            zmin = temp;
        if ((temp = addr->edges[nedge]->v2->z) > zmax)
            zmax = temp;
        if (temp < zmin)
            zmin = temp;
    }

    addr->xmin = xmin;
    addr->xmax = xmax;
    addr->ymin = ymin;
    addr->ymax = ymax;
    addr->zmin = zmin;
    addr->zmax = zmax;

    printf("POLY LIMITS p%d  x %f %f  y %f %f  z %f %f\n", addr->id_poly,
            xmin, xmax, ymin, ymax, zmin, zmax);

    return (0);

}

/*** function to set solid limits */

int set_solid_limits(struct solid *addr) {
    int npoly;
    double xmin, ymin, zmin, xmax, ymax, zmax, temp;

    xmin = ymin = zmin = REAL_MAX;
    xmax = ymax = zmax = -REAL_MAX;

    /* check limits for each polygon */

    for (npoly = 0; npoly < addr->n_poly; npoly++) {
        if ((temp = addr->poly[npoly]->xmax) > xmax)
            xmax = temp;
        if ((temp = addr->poly[npoly]->xmin) < xmin)
            xmin = temp;
        if ((temp = addr->poly[npoly]->ymax) > ymax)
            ymax = temp;
        if ((temp = addr->poly[npoly]->ymin) < ymin)
            ymin = temp;
        if ((temp = addr->poly[npoly]->zmax) > zmax)
            zmax = temp;
        if ((temp = addr->poly[npoly]->zmin) < zmin)
            zmin = temp;
    }

    addr->xmin = xmin;
    addr->xmax = xmax;
    addr->ymin = ymin;
    addr->ymax = ymax;
    addr->zmin = zmin;
    addr->zmax = zmax;

    printf("SOLID LIMITS s%d  x %f %f  y %f %f  z %f %f\n", addr->id_solid,
            xmin, xmax, ymin, ymax, zmin, zmax);

    return (0);

}

/*** function to calculate velocity from  polygon model for 2D to 3D case */

double get_poly_vel_2D3D(double xpos, double ypos, double zpos, char wavetype,
        double* pdensity, int iden, int *pimodel) {
    double xtrans, ytrans;


    /* 3D to 2D conversion */

    xtrans = xpos - Mod2D3D_origx;
    ytrans = ypos - Mod2D3D_origy;

    xtrans = xtrans * Mod2D3D_cosang - ytrans * Mod2D3D_sinang;

    return (get_poly_vel(xtrans, zpos, wavetype, pdensity, iden, pimodel));

}

/*** function to calculate velocity from polygon model */

double get_poly_vel(double xpos, double zpos, char wavetype, double* density,
        int iden, int *imodel) {
    double vel;
    struct polygon *addr;


    if ((addr = poly_head) == NULL) {
        return (-1.0);
    }



    do {
        if (xpos >= addr->xmin && xpos <= addr->xmax &&
                zpos >= addr->zmin && zpos <= addr->zmax) {

            if (inside_poly(xpos, zpos, addr)) {

                if (wavetype == 'P')
                    vel = addr->vpref +
                        addr->vpgrad * (zpos - addr->ref);
                else
                    vel = addr->vsref +
                        addr->vsgrad * (zpos - addr->ref);

                if (iden == 1)
                    *density = addr->denref +
                        addr->dengrad * (zpos - addr->ref);

                /* set model element code */

                *imodel = POLYOFFSET + addr->id_poly;

                return (vel);
            }
        }
    } while ((addr = addr->next) != poly_head);

    return (-1.0);
}

/*** function to determine if a point is inside a simple polygon */

int inside_poly(double xpt, double zpt, struct polygon *addr)
/*!double xpt, zpt;
struct polygon *addr;*/ {
    int nedge;

    for (nedge = 0; nedge < addr->n_edges; nedge++) {
        /*		printf("TEST e%d x %f y %f test=%f\n", addr->edges[nedge]->id_edge,
                                xpt, zpt,
                                xpt * addr->norm[nedge].a + zpt * addr->norm[nedge].b +
                                addr->norm[nedge].c);
         */
        if (xpt * addr->norm[nedge].a + zpt * addr->norm[nedge].b +
                addr->norm[nedge].c < 0.0)
            return (0);
    }

    return (1);
}

/*** function to calculate velocity from polygon model */

double get_solid_vel(double xpos, double ypos, double zpos, char wavetype, double *density, int iden) {
    double vel;
    struct solid *addr;


    if ((addr = solid_head) == NULL) {
        return (-1.0);
    }

    do {
        if (xpos >= addr->xmin && xpos <= addr->xmax &&
                ypos >= addr->ymin && ypos <= addr->ymax &&
                zpos >= addr->zmin && zpos <= addr->zmax) {
            if (inside_solid(xpos, ypos, zpos, addr)) {
                if (wavetype == 'P')
                    vel = addr->vpref +
                        addr->vpgrad * (zpos - addr->ref);
                else
                    vel = addr->vsref +
                        addr->vsgrad * (zpos - addr->ref);
                if (iden == 1)
                    *density = addr->denref +
                        addr->dengrad * (zpos - addr->ref);
                return (vel);
            }
        }
    } while ((addr = addr->next) != solid_head);

    return (-1.0);
}

/*** function to determine if a point is inside a simple solid */

int inside_solid(double xpt, double ypt, double zpt, struct solid *addr)
/*!double xpt, ypt, zpt;
struct solid *addr;*/ {
    int npoly;

    for (npoly = 0; npoly < addr->n_poly; npoly++) {
        /*		printf("TEST p%d x %f y %f z %f test=%f\n", addr->poly[npoly]->id_poly,
                                xpt, ypt, zpt,
                                xpt * addr->norm[npoly].a + ypt * addr->norm[npoly].b +
                                zpt * addr->norm[npoly].c +
                                addr->norm[npoly].d);*/

        if (xpt * addr->norm[npoly].a + ypt * addr->norm[npoly].b +
                zpt * addr->norm[npoly].c + addr->norm[npoly].d < 0.0)
            return (0);
    }

    return (1);
}

/* function to read finite difference grid model */

void read_fdiff_vel(char* fname) {
    FILE *fp_grid;
    int gridsize, nz, nx;
    double vmean, smean;
    float *addr;

    if ((fp_grid = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open velocity grid file:\n");
        fprintf(stderr, "  %s\n", fname);
    }

    fseek(fp_grid, (long) sizeof (int), SEEK_SET);
    fread(&fdgrid_numx, sizeof (int), 1, fp_grid);
    fread(&fdgrid_numz, sizeof (int), 1, fp_grid);
    fdgrid_numx++;
    fdgrid_numz++;
    printf("Finite Diff Vel grid:\n  Nx %d  Nz %d\n",
            fdgrid_numx, fdgrid_numz);
    fdgrid_xstep = (fdgrid_xmax - fdgrid_xmin) / (double) (fdgrid_numx - 1);
    fdgrid_zstep = (fdgrid_zmax - fdgrid_zmin) / (double) (fdgrid_numz - 1);

    gridsize = fdgrid_numx * fdgrid_numz * sizeof (double);
    if ((fdgrid_array = (float *) malloc((size_t) gridsize)) == NULL)
        fprintf(stderr,
            "ERROR: Cannot allocate array for grid velocities.\n");

    fseek(fp_grid, (long) (2 * sizeof (int)), SEEK_CUR);
    for (nz = 0; nz < fdgrid_numz; nz++) {
        fread(fdgrid_array + (nz * fdgrid_numx), sizeof (float),
                fdgrid_numx, fp_grid);
        if (nz == 0 || nz == (fdgrid_numz - 1))
            printf("  Row nz = %4d: %f  %f  ...  %f  %f\n", nz, *(fdgrid_array + (nz * fdgrid_numx)), *(fdgrid_array + (nz * fdgrid_numx) + 1), *(fdgrid_array + (nz * fdgrid_numx) + fdgrid_numx - 2), *(fdgrid_array + (nz * fdgrid_numx) + fdgrid_numx - 1));
    }

    vmean = (smean = 0.0);
    for (nz = 0; nz < fdgrid_numz; nz++) {
        for (nx = 0; nx < fdgrid_numx; nx++) {
            addr = fdgrid_array + (nz * fdgrid_numx + nx);
            *addr /= 1000.0f;
            vmean += *addr;
            *addr = 1.0f / *addr;
            smean += *addr;
        }
    }
    vmean /= fdgrid_numx * fdgrid_numz;
    smean /= fdgrid_numx * fdgrid_numz;
    printf("  vmean %lf  1/smean %lf\n", vmean, 1.0 / smean);
    vmodel_vmean = vmean;

}

/* function to calculate velocity from finite difference grid model */

double get_fdiff_vel(double xpos, double zpos, char wavetype, double *density, int iden) {
    int zindex, xindex;
    double vel, slow;

    /*	if (zpos < fdgrid_zmin)
                    zpos = fdgrid_zmin + (fdgrid_zmin - zpos);
            if (xpos < fdgrid_xmin)
                    xpos = fdgrid_xmin + (fdgrid_xmin - xpos);
            if (zpos > fdgrid_zmax)
                    zpos = fdgrid_zmax - (zpos - fdgrid_zmax);
            if (xpos > fdgrid_xmax)
                    xpos = fdgrid_xmax - (xpos - fdgrid_xmax);
     */
    if (zpos < fdgrid_zmin)
        zpos = fdgrid_zmin;
    if (xpos < fdgrid_xmin)
        xpos = fdgrid_xmin;
    if (zpos > fdgrid_zmax)
        zpos = fdgrid_zmax;
    if (xpos > fdgrid_xmax)
        xpos = fdgrid_xmax;


    if (zpos >= fdgrid_zmin && zpos <= fdgrid_zmax &&
            xpos >= fdgrid_xmin && xpos <= fdgrid_xmax) {
        zindex = (int) ((zpos - fdgrid_zmin) / fdgrid_zstep);
        xindex = (int) ((xpos - fdgrid_xmin) / fdgrid_xstep);
        slow = *(fdgrid_array + (zindex * fdgrid_numx + xindex));
        vel = 1.0 / slow;
    } else
        vel = vmodel_vmean;

    if (iden)
        *density = 2.7;
    return (vel);
}


