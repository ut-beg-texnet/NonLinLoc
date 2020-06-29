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


// function to run wavefront-ray algorithm green3d using a function call

int RunGreen3d(GridDesc* pmgrid, SourceDesc* psource, GridDesc* ptt_grid)
{

	// Wavefront ray tracing (green3d) declarations
	double km2m = 1000.0;
	float dxm, dym, dzm;
	float xlimits[2], ylimits[2], zlimits[2], xsrc[3];
	int ntab;
	float *tab;
	int NFRONT, ND , NK, NPAR;


		// load parameters

		dxm = km2m * pmgrid->dx;
		dym = km2m * pmgrid->dy;
		dzm = km2m * pmgrid->dz;
		xlimits[0] = km2m * pmgrid->origx + km2m * pmgrid->dx;
		xlimits[1] = km2m * pmgrid->origx + 
			(double) (pmgrid->numx - 2) * km2m * pmgrid->dx;
		ylimits[0] = km2m * pmgrid->origy + km2m * pmgrid->dy;
		ylimits[1] = km2m * pmgrid->origy + 
			(double) (pmgrid->numy - 2) * km2m * pmgrid->dy;
		zlimits[0] = km2m * pmgrid->origz + km2m * pmgrid->dz;
		zlimits[1] = km2m * pmgrid->origz + 
			(double) (pmgrid->numz - 2) * km2m * pmgrid->dz;

		wvfrnt_targ_orient[0][0] = 
				km2m * (ptt_grid->origy + ptt_grid->dy);
		wvfrnt_targ_orient[0][1] = 
				km2m * (ptt_grid->origx + ptt_grid->dx);
		wvfrnt_targ_orient[0][2] = 
				km2m * (ptt_grid->origz + ptt_grid->dz);
		wvfrnt_targ_orient[1][0] = km2m * ptt_grid->dy;
		wvfrnt_targ_orient[1][1] = 0.0;
		wvfrnt_targ_orient[1][2] = 0.0;
		wvfrnt_targ_orient[2][0] = 0.0;
		wvfrnt_targ_orient[2][1] = km2m * ptt_grid->dx;
		wvfrnt_targ_orient[2][2] = 0.0;
		wvfrnt_targ_orient[3][0] = 0.0;
		wvfrnt_targ_orient[3][1] = 0.0;
		wvfrnt_targ_orient[3][2] = km2m * ptt_grid->dz;
      		// call inv3x3(targ_orient(1,2),targ_orient(1,5)) [modelrese.f]
		(void) inv3x3_(wvfrnt_targ_orient[1], wvfrnt_targ_orient[4]);

		xsrc[0] = km2m * psource->y;
		xsrc[1] = km2m * psource->x;
		xsrc[2] = km2m * psource->z;

		// tab memory allocation
		// parameter (nfront=2,nd=2,nk=2,npar=48)
		NFRONT = 2; ND = 2; NK = 2; NPAR = 48;
		ntab = 1000 * (NFRONT*(NPAR+NK+2*ND+30));
		if ((tab = (float *) malloc((size_t) (ntab * sizeof(float)))) 
				== NULL) {
			nll_puterr(
"ERROR allocating memory for wavefront/green3a \"tab\" buffer.\n");
			exit(EXIT_ERROR_MEMORY);
		}

 		// call wavefront algorithm
		// Note: green3a map storage is FORTRAN (z,x,y) = C (y,x,z)
		//   but we use C (x,y,z), thus we must exchange some
		//   x and y arguments to green3a.
     		(void) green3a_(
			&(pmgrid->numy), &(pmgrid->numx), &(pmgrid->numz),
			&dym, &dxm, &dzm,
			ylimits, xlimits, zlimits, pmgrid->buffer,
 			&(ptt_grid->numy), &(ptt_grid->numx), &(ptt_grid->numz),
			&wvfrnt_nir, &wvfrnt_npr, wvfrnt_targ_orient,
			ptt_grid->buffer, wvfrnt_imap, xsrc,
			&wvfrnt_fi1min, &wvfrnt_fi2min, 
			&wvfrnt_fi1max, &wvfrnt_fi2max,
			&wvfrnt_dxmin2, &wvfrnt_dpmin2, &wvfrnt_dtemps,
			&ntab, tab);


	return(0);


}


