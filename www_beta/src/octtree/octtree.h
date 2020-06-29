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


/*  octree.h

	include file for octree search

*/



/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

	/* misc defines */

#ifndef SMALL_DOUBLE
#define SMALL_DOUBLE 1.0e-20
#endif
#ifndef LARGE_DOUBLE
#define LARGE_DOUBLE 1.0e20
#endif
#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif
#ifndef VERY_LARGE_DOUBLE
#define VERY_LARGE_DOUBLE 1.0e30
#endif


/*------------------------------------------------------------/ */
/* structures */
/*------------------------------------------------------------/ */

#define VALUE_IS_LOG_PROB_DENSITY_IN_NODE 0     // log of pdf (does not take into account volume of cell)
#define VALUE_IS_PROB_DENSITY_IN_NODE 1     // pdf (does not take into account volume of cell)
#define VALUE_IS_PROBABILITY_IN_NODE 2     // pdf * cell volume (takes into account volume of cell)

/* octree node */

typedef struct octnode* OctNodePtr;
typedef struct octnode
{
	OctNodePtr parent;		/* parent node */
	Vect3D center;			/* absolute coordinates of center */
	Vect3D ds;			/* length of sides */
	int level;                      // level of node in oect-tree hierarchy (0 = top, largest)
	double value;			/* node value */
	OctNodePtr child[2][2][2];	/* child nodes */
	char isLeaf;			/* leaf flag, 1=leaf */
	void *pdata;		/* additional data */
} OctNode;



/* 3D tree with Nx, Ny, Nz arbitrary */

typedef struct
{
	OctNode**** nodeArray;      // parent nodes
	int data_code;              // data type code, application dependent
	int numx, numy, numz;       // grid size
 	Vect3D orig;                // orig (km)
	Vect3D ds;                  // len side (km) (nominal for ds.x for spherical case)
        double* ds_x;               // array of true ds.x values for spherical case
        int* num_x;                 // array of true num_x values for spherical case
	double integral;
        int isSpherical;            // =1 if Tree3D is spherical, 0 otherwise
}
Tree3D;


/* structure for storing results */

typedef struct resultTreeNode* ResultTreeNodePtr;
typedef struct resultTreeNode {
	ResultTreeNodePtr left;		/* address of left node */
	ResultTreeNodePtr right;	/* address of right node */
	double value;			/* sort value */
	int level;			/* level of node in oect-tree hierarchy (0 = top, largest) */
	double volume;		/* volume, node volume depends on geometry in physical space, may not be dx*dy*dz */
	OctNode* pnode;			/* corresponding octree node */
} ResultTreeNode;



/* */
/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* globals  */
/*------------------------------------------------------------/ */

//EXTERN_TXT char fn_control[MAXLINE];	/* control file name */

/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

Tree3D* newTree3D(int data_code, int numx, int numy, int numz,
	double origx, double origy, double origz,
	double dx,  double dy,  double dz, double value, double integral, void *pdata);
Tree3D* newTree3D_spherical(int data_code, int numx_nominal, int numy, int numz,
        double origx, double origy, double origz,
        double dx_nominal, double dy, double dz, double value, double integral, void *pdata);
double get_dx_spherical(double dx_nominal, double origx, double x_max, double center_y, int *pnum_x);
OctNode* newOctNode(OctNode* parent, Vect3D center, Vect3D ds, double value, void *pdata);
void subdivide(OctNode* parent, double value, void *pdata);
void freeTree3D(Tree3D* tree, int freeDataPointer);
void freeNode(OctNode* node, int freeDataPointer);
OctNode* getLeafNodeContaining(Tree3D* tree, Vect3D coords);
OctNode* getLeafContaining(OctNode* node, double x, double y, double z);

ResultTreeNode* addResult(ResultTreeNode* prtn, double value, double volume, OctNode* pnode);
void freeResultTree(ResultTreeNode* prtn);
ResultTreeNode* getHighestValue(ResultTreeNode* prtn);
ResultTreeNode* getHighestLeafValue(ResultTreeNode* prtree);
ResultTreeNode* getHighestLeafValueMinSize(ResultTreeNode* prtree, double sizeMinX, double sizeMinY, double sizeMinZ);
ResultTreeNode* getHighestLeafValueLESpecifiedSize(ResultTreeNode* prtree, double sizeX, double sizeY, double sizeZ);
ResultTreeNode* getHighestLeafValueOfSpecifiedSize(ResultTreeNode* prtree, double sizeX, double sizeY, double sizeZ);
ResultTreeNode* getHighestLeafValueAtSpecifiedLevel(ResultTreeNode* prtree, int level);

Tree3D* readTree3D(FILE *fpio);
int readNode(FILE *fpio, OctNode* node);
int writeTree3D(FILE *fpio, Tree3D* tree);
int writeNode(FILE *fpio, OctNode* node);

int nodeContains(OctNode* node, double x, double y, double z);
int extendedNodeContains(OctNode* node, double x, double y, double z, int checkZ);

int getScatterSampleResultTreeAtLevels(ResultTreeNode* prtree, int value_type, int num_scatter,
        double integral, float* fdata, int npoints, int* pfdata_index,
        double oct_node_value_ref, double *poct_tree_scatter_volume, int level_min, int level_max);
int getScatterSampleResultTree(ResultTreeNode* prtree, int value_type, int num_scatter,
        double integral, float* fdata, int npoints, int* pfdata_index,
        double oct_node_value_max, double *poct_tree_scatter_volume);
double convertOcttreeValuesToProb(ResultTreeNode* prtree, double sum, double oct_node_value_max);
double integrateResultTreeAtLevels(ResultTreeNode* prtree, int value_type, double sum, double oct_node_value_max, int level_min, int level_max);
double integrateResultTree(ResultTreeNode* prtree, int value_type, double sum, double oct_node_value_max);
ResultTreeNode* createResultTree(ResultTreeNode* prtree, ResultTreeNode* pnew_rtree);


/* */
/*------------------------------------------------------------/ */


