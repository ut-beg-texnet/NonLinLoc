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


/*   octree.c

        octree search functions

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    02NOV2000  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#define EXTERN_MODE 1

#include "geometry.h"
#include "ran1.h"
#include "octtree.h"

/*** function to create a new OctNode */

OctNode* newOctNode(OctNode* parent, Vect3D center, Vect3D ds, double value, void *pdata) {

    int l, m, n;
    OctNode* node;

    node = (OctNode*) malloc(sizeof (OctNode));

    node->parent = parent;
    node->center = center;
    node->ds = ds;
    if (parent != NULL)
        node->level = parent->level + 1;
    else
        node->level = 0;
    node->value = value;
    node->pdata = pdata;

    for (l = 0; l < 2; l++)
        for (m = 0; m < 2; m++)
            for (n = 0; n < 2; n++)
                node->child[l][m][n] = NULL;

    node->isLeaf = 1;

    return (node);
}

/*** function to free a Tree3D ***/

void freeTree3D(Tree3D* tree, int freeDataPointer) {

    // AEH/AJL 20080709
    if (tree == NULL) return;

    if (tree->isSpherical) {
        free(tree->ds_x);
        free(tree->num_x);
    }

    int ix, iy, iz;

    for (ix = 0; ix < tree->numx; ix++) {
        for (iy = 0; iy < tree->numy; iy++) {
            for (iz = 0; iz < tree->numz; iz++) {
                if (tree->nodeArray[ix][iy][iz] != NULL) // case of Tree3D_spherical
                    freeNode(tree->nodeArray[ix][iy][iz], freeDataPointer);
            }
            free(tree->nodeArray[ix][iy]);
        }
        free(tree->nodeArray[ix]);
    }

    // AJL - 20080710  (valgrind)
    free(tree->nodeArray);

    free(tree);

}

/*** function to create a new Tree3D - an x, y, z array of octtree root nodes ***/

Tree3D* newTree3D(int data_code, int numx, int numy, int numz,
        double origx, double origy, double origz,
        double dx, double dy, double dz, double value, double integral, void *pdata) {

    int ix, iy, iz;
    OctNode**** garray;
    Tree3D* tree;
    Vect3D center, ds;

    // create Tree3D structure
    tree = (Tree3D*) malloc(sizeof (Tree3D));
    if (tree == NULL) {
        return (NULL);
    }
    // alocate node array in x, y and z
    garray = (OctNode ****) malloc((size_t) numx * sizeof (OctNode***));
    if (garray == NULL) {
        free(tree);
        return (NULL);
    }
    tree->ds_x = NULL;
    tree->num_x = NULL;

    ds.x = dx;
    ds.y = dy;
    ds.z = dz;

    for (ix = 0; ix < numx; ix++) {
        center.x = origx + (double) ix * dx + dx / 2.0;
        if ((garray[ix] = (OctNode ***) malloc((size_t) numy *
                sizeof (OctNode**))) == NULL)
            return (NULL);
        for (iy = 0; iy < numy; iy++) {
            center.y = origy + (double) iy * dy + dy / 2.0;
            if ((garray[ix][iy] = (OctNode **) malloc((size_t) numz *
                    sizeof (OctNode*))) == NULL)
                return (NULL);
            for (iz = 0; iz < numz; iz++) {
                center.z = origz + (double) iz * dz + dz / 2.0;
                garray[ix][iy][iz] = newOctNode(NULL, center, ds, value, pdata);
            }
        }
    }


    // set Tree3D structure values
    tree->nodeArray = garray;
    tree->data_code = data_code;
    tree->numx = numx;
    tree->numy = numy;
    tree->numz = numz;
    tree->orig.x = origx;
    tree->orig.y = origy;
    tree->orig.z = origz;
    tree->ds = ds;
    tree->integral = integral;
    tree->isSpherical = 0;

    return (tree);

}

#define DE2RA (M_PI/180.0)

/*** function to create a new Tree3D - an x, y, z array of octtree root nodes
 * creates a Tree3D for a sphere:
 *    width of cells in azimuth ior longitude is approximately constant,
 *    giving fewer azimuth cells with increasing inclination or latitiude (y)
 ***/

Tree3D* newTree3D_spherical(int data_code, int numx_nominal, int numy, int numz,
        double origx, double origy, double origz,
        double dx_nominal, double dy, double dz, double value, double integral, void *pdata) {

    int ix, iy, iz;
    OctNode**** garray;
    Tree3D* tree;
    Vect3D center, ds;


    // create Tree3D structure
    tree = (Tree3D*) malloc(sizeof (Tree3D));
    if (tree == NULL) {
        return (NULL);
    }
    // alocate node array in x, y and z
    garray = (OctNode ****) malloc((size_t) numx_nominal * sizeof (OctNode***));
    if (garray == NULL) {
        free(tree);
        return (NULL);
    }
    // alocate ds_x array
    tree->ds_x = (double*) malloc((size_t) numy * sizeof (double));
    if (tree->ds_x == NULL) {
        free(garray);
        free(tree);
        return (NULL);
    }
    // alocate num_x array
    tree->num_x = (int*) malloc((size_t) numy * sizeof (int));
    if (tree->num_x == NULL) {
        free(garray);
        free(tree->ds_x);
        free(tree);
        return (NULL);
    }

    ds.x = dx_nominal;
    ds.y = dy;
    ds.z = dz;
    Vect3D tree_ds = ds;

    double x_max = origx + (double) numx_nominal * dx_nominal;
    double dx;
    int num_x;

    for (ix = 0; ix < numx_nominal; ix++) {
        if ((garray[ix] = (OctNode ***) malloc((size_t) numy *
                sizeof (OctNode**))) == NULL)
            return (NULL);
        for (iy = 0; iy < numy; iy++) {
            center.y = origy + (double) iy * dy + dy / 2.0;
            if ((garray[ix][iy] = (OctNode **) malloc((size_t) numz *
                    sizeof (OctNode*))) == NULL)
                return (NULL);
            dx = get_dx_spherical(dx_nominal, origx, x_max, center.y, &num_x);
            //printf("dx %g  dx_nominal %g  origx %g  x_max %g  center.y %g  num_x %d\n", dx, dx_nominal, origx, x_max, center.y, num_x);
            if (ix == 0) {
                tree->ds_x[iy] = dx;
                tree->num_x[iy] = num_x;
            }
            for (iz = 0; iz < numz; iz++) {
                if (ix < num_x) {
                    ds.x = dx;
                    center.x = origx + (double) ix * dx + dx / 2.0;
                    center.z = origz + (double) iz * dz + dz / 2.0;
                    garray[ix][iy][iz] = newOctNode(NULL, center, ds, value, pdata);
                } else {
                    garray[ix][iy][iz] = NULL;
                }
            }
        }
    }


    // set Tree3D structure values
    tree->nodeArray = garray;
    tree->data_code = data_code;
    tree->numx = numx_nominal;
    tree->numy = numy;
    tree->numz = numz;
    tree->orig.x = origx;
    tree->orig.y = origy;
    tree->orig.z = origz;
    tree->ds = tree_ds;
    tree->integral = integral;
    tree->isSpherical = 1;

    return (tree);

}

/*** funciton to calculate dx values as a funciton of y value for a Tree3D spherical*/

double get_dx_spherical(double dx_nominal, double origx, double x_max, double center_y, int *pnum_x) {

    double dx = dx_nominal / cos(center_y * DE2RA);
    *pnum_x = (int) (0.5 + (x_max - origx) / dx);
    //if (*pnum_x < 4)
    //    *pnum_x = 4;
    dx = (x_max - origx) / (double) (*pnum_x); // correct x step for integer rounding
    return (dx);

}

/*** function to sudivide a node into child nodes ***/

void subdivide(OctNode* parent, double value, void *pdata) {

    int ix, iy, iz;
    Vect3D center, ds;

    ds.x = parent->ds.x / 2.0;
    ds.y = parent->ds.y / 2.0;
    ds.z = parent->ds.z / 2.0;

    // alocate nodes in x, y and z

    for (ix = 0; ix < 2; ix++) {
        center.x = parent->center.x + (double) (2 * ix - 1) * ds.x / 2.0;
        for (iy = 0; iy < 2; iy++) {
            center.y = parent->center.y + (double) (2 * iy - 1) * ds.y / 2.0;
            for (iz = 0; iz < 2; iz++) {
                center.z = parent->center.z + (double) (2 * iz - 1) * ds.z / 2.0;
                parent->child[ix][iy][iz] = newOctNode(parent, center, ds, value, pdata);
            }
        }
    }

    if (parent != NULL)
        parent->isLeaf = 0;

}

/*** function to free an OctNode and all its child nodes ***/

void freeNode(OctNode* node, int freeDataPointer) {

    int ix, iy, iz;
    for (ix = 0; ix < 2; ix++) {
        for (iy = 0; iy < 2; iy++) {
            for (iz = 0; iz < 2; iz++) {
                if (node->child[ix][iy][iz] != NULL)
                    freeNode(node->child[ix][iy][iz], freeDataPointer);
            }
        }
    }

    // try to free data
    if (freeDataPointer)
        free(node->pdata);
    free(node);

}

/*** function to get the leaf node in a Tree3D containing the given x, y, z coordinates ***/

OctNode* getLeafNodeContaining(Tree3D* tree, Vect3D coords) {

    int ix, iy, iz;
    OctNode* rootNode;

    // get indices of root in tree
    // z
    iz = (int) ((coords.z - tree->orig.z) / tree->ds.z);
    if (iz < 0 || iz >= tree->numz)
        return (NULL);
    // y
    iy = (int) ((coords.y - tree->orig.y) / tree->ds.y);
    if (iy < 0 || iy >= tree->numy)
        return (NULL);
    // x
    if (tree->isSpherical) {
        double dx = tree->ds_x[iy];
        int num_x = tree->num_x[iy];
        ix = (int) ((coords.x - tree->orig.x) / dx);
        if (ix < 0 || ix >= num_x)
            return (NULL);
    } else {
        ix = (int) ((coords.x - tree->orig.x) / tree->ds.x);
        if (ix < 0 || ix >= tree->numx)
            return (NULL);
    }

    rootNode = tree->nodeArray[ix][iy][iz];

    if (rootNode == NULL)
        return (NULL);

    return (getLeafContaining(rootNode, coords.x, coords.y, coords.z));
}

/*** function to get the leaf node in a node containing the given x, y, z coordinates ***/

OctNode* getLeafContaining(OctNode* node, double x, double y, double z) {

    int ix = 0, iy = 0, iz = 0;
    OctNode* childNode;

    // get indices of child in this node
    if (x >= node->center.x)
        ix = 1;
    if (y >= node->center.y)
        iy = 1;
    if (z >= node->center.z)
        iz = 1;

    childNode = node->child[ix][iy][iz];

    if (childNode == NULL)
        return (node);

    return (getLeafContaining(childNode, x, y, z));

}

/*** function to check if node contains the given x, y, z coordinates ***/

int nodeContains(OctNode* node, double x, double y, double z) {

    Vect3D ds;
    double dx2, dy2, dz2;

    ds = node->ds;

    dx2 = ds.x / 2.0;
    if (x < node->center.x - dx2 || x > node->center.x + dx2)
        return (0);

    dy2 = ds.y / 2.0;
    if (y < node->center.y - dy2 || y > node->center.y + dy2)
        return (0);

    dz2 = ds.z / 2.0;
    if (z < node->center.z - dz2 || z > node->center.z + dz2)
        return (0);

    return (1);

}

/*** function to check if node or half-segments of adjacent equal-size nodes contains the given x, y, z coordinates ***/

int extendedNodeContains(OctNode* node, double x, double y, double z, int checkZ) {

    Vect3D ds;
    double dx2, dy2, dz2;

    ds = node->ds;

    dx2 = ds.x;
    if (x < node->center.x - dx2 || x > node->center.x + dx2)
        return (0);

    dy2 = ds.y;
    if (y < node->center.y - dy2 || y > node->center.y + dy2)
        return (0);

    if (checkZ) {
        dz2 = ds.z;
        if (z < node->center.z - dz2 || z > node->center.z + dz2)
            return (0);
    }

    return (1);

}

/*** function to put Octtree node in results tree in order of value */

ResultTreeNode* addResult(ResultTreeNode* prtree, double value, double volume, OctNode* pnode) {

    // put address in result tree based on value

    if (prtree == NULL) { // at empty node
        if ((prtree = (ResultTreeNode*) malloc(sizeof (ResultTreeNode))) == NULL)
            fprintf(stderr, "ERROR allocating memory for result-tree node.\n");
        prtree->value = value;
        prtree->level = pnode->level;
        prtree->volume = volume; // node volume depends on geometry in physical space, may not be dx*dy*dz
        prtree->pnode = pnode;
        prtree->left = prtree->right = NULL;

        /* 20101213 AJL removed this block - causes memory errors in trace_processing/location.c
           } else if (value == prtree->value) { // prevent assymetric tree if multiple identical values
            if (get_rand_int(-10000, 9999) < 0)
                prtree->left = addResult(prtree->left, value, volume, pnode);
            else
                prtree->right = addResult(prtree->right, value, volume, pnode);
         */

    } else if (value < prtree->value) {
        prtree->left = addResult(prtree->left, value, volume, pnode);

    } else {
        prtree->right = addResult(prtree->right, value, volume, pnode);
    }

    return (prtree);
}

/*** function to free results tree */

void freeResultTree(ResultTreeNode* prtree) {

    if (prtree == NULL)
        return;

    if (prtree->left != NULL)
        freeResultTree(prtree->left);
    if (prtree->right != NULL)
        freeResultTree(prtree->right);
    free(prtree);
}

/*** function to get ResultTreeNode with highest value */

ResultTreeNode* getHighestValue(ResultTreeNode* prtree) {

    if (prtree == NULL)
        return (NULL);

    if (prtree->right == NULL) /* right child is empty */
        return (prtree);

    return (getHighestValue(prtree->right));
}

/*** function to get ResultTree Leaf Node with highest value */

ResultTreeNode* getHighestLeafValue(ResultTreeNode* prtree) {

    if (prtree == NULL)
        return (NULL);

    ResultTreeNode* prtree_returned = NULL;

    if (prtree->right != NULL)
        prtree_returned = getHighestLeafValue(prtree->right);

    if (prtree_returned != NULL) // right leaf descendent
        return (prtree_returned); // thus highest value leaf
    else // no right leaf descendents
        if (prtree->pnode->isLeaf) // this is leaf
        return (prtree); // thus highest value leaf

    if (prtree->left != NULL) // look for left leaf descendents
        prtree_returned = getHighestLeafValue(prtree->left);

    return (prtree_returned);
}

/*** function to get ResultTree Leaf Node with highest value, node must be larger than specified minimum size */

ResultTreeNode* getHighestLeafValueMinSize(ResultTreeNode* prtree, double sizeMinX, double sizeMinY, double sizeMinZ) {

    if (prtree == NULL) { // no nodes in tree
        return (NULL);
    }

    ResultTreeNode* prtree_returned = NULL;

    if (prtree->right != NULL)
        prtree_returned = getHighestLeafValueMinSize(prtree->right, sizeMinX, sizeMinY, sizeMinZ);

    if (prtree_returned != NULL) // right leaf descendent
        return (prtree_returned); // thus highest value leaf
    else { // no right leaf descendents
        if (prtree->pnode->isLeaf // this is leaf
                && prtree->pnode->ds.x >= sizeMinX
                && prtree->pnode->ds.y >= sizeMinY
                && prtree->pnode->ds.z >= sizeMinZ) // not too small
            return (prtree); // thus highest value leaf that is not too small
    }

    if (prtree->left != NULL) // look for left leaf descendents
        prtree_returned = getHighestLeafValueMinSize(prtree->left, sizeMinX, sizeMinY, sizeMinZ);

    return (prtree_returned);
}


#define SIZE_TOLERANCE 1.0e-20

/*** function to get ResultTree Leaf Node with highest value, node size must be less than or equal to specified maximum size */

ResultTreeNode* getHighestLeafValueLESpecifiedSize(ResultTreeNode* prtree, double sizeX, double sizeY, double sizeZ) {

    ResultTreeNode* prtree_returned = NULL;

    if (prtree->right != NULL)
        prtree_returned = getHighestLeafValueLESpecifiedSize(prtree->right, sizeX, sizeY, sizeZ);

    if (prtree_returned != NULL) // right leaf descendent
        return (prtree_returned); // thus highest value leaf
    else { // no right leaf descendents
        if (prtree->pnode->isLeaf // this is leaf
                && ((sizeX < 0.0) || (prtree->pnode->ds.x - sizeX < SIZE_TOLERANCE))
                && ((sizeY < 0.0) || (prtree->pnode->ds.y - sizeY < SIZE_TOLERANCE))
                && ((sizeZ < 0.0) || (prtree->pnode->ds.z - sizeZ < SIZE_TOLERANCE))
                ) // is specified size
            return (prtree); // thus highest value leaf that is not too small
    }

    if (prtree->left != NULL) // look for left leaf descendents
        prtree_returned = getHighestLeafValueLESpecifiedSize(prtree->left, sizeX, sizeY, sizeZ);

    return (prtree_returned);
}

/*** function to get ResultTree Leaf Node with highest value, node size must be equal to specified maximum size */

ResultTreeNode* getHighestLeafValueOfSpecifiedSize(ResultTreeNode* prtree, double sizeX, double sizeY, double sizeZ) {

    ResultTreeNode* prtree_returned = NULL;

    if (prtree->right != NULL)
        prtree_returned = getHighestLeafValueOfSpecifiedSize(prtree->right, sizeX, sizeY, sizeZ);

    if (prtree_returned != NULL) // right leaf descendent
        return (prtree_returned); // thus highest value leaf
    else { // no right leaf descendents
        if (prtree->pnode->isLeaf // this is leaf
                && ((sizeX < 0.0) || (fabs(prtree->pnode->ds.x - sizeX) < SIZE_TOLERANCE))
                && ((sizeY < 0.0) || (fabs(prtree->pnode->ds.y - sizeY) < SIZE_TOLERANCE))
                && ((sizeZ < 0.0) || (fabs(prtree->pnode->ds.z - sizeZ) < SIZE_TOLERANCE))
                ) // is specified size
            return (prtree); // thus highest value leaf that is not too small
    }

    if (prtree->left != NULL) // look for left leaf descendents
        prtree_returned = getHighestLeafValueOfSpecifiedSize(prtree->left, sizeX, sizeY, sizeZ);

    return (prtree_returned);
}

/*** function to get ResultTree Leaf Node at specified level with highest value */

ResultTreeNode* getHighestLeafValueAtSpecifiedLevel(ResultTreeNode* prtree, int level) {

    ResultTreeNode* prtree_returned = NULL;

    if (prtree->right != NULL)
        prtree_returned = getHighestLeafValueAtSpecifiedLevel(prtree->right, level);

    if (prtree_returned != NULL) // right leaf descendent
        return (prtree_returned); // thus highest value leaf
    else { // no right leaf descendents
        if ((prtree->level == level) // is specified level
                && prtree->pnode->isLeaf) // this is leaf
            return (prtree); // thus highest value leaf that is at specfied level
    }

    if (prtree->left != NULL) // look for left leaf descendents
        prtree_returned = getHighestLeafValueAtSpecifiedLevel(prtree->left, level);

    return (prtree_returned);
}

/*** function to read a Tree3D ***/

// TODO: modify to support spherical Tree3D!

Tree3D* readTree3D(FILE *fpio) {

    int istat;
    int istat_cum;
    int ix, iy, iz;

    Tree3D* tree = NULL;
    int data_code;
    int numx, numy, numz;
    Vect3D orig;
    Vect3D ds;
    double integral;


    istat = fread(&data_code, sizeof (int), 1, fpio);
    istat += fread(&numx, sizeof (int), 1, fpio);
    istat += fread(&numy, sizeof (int), 1, fpio);
    istat += fread(&numz, sizeof (int), 1, fpio);
    istat += fread(&orig, sizeof (Vect3D), 1, fpio);
    istat += fread(&ds, sizeof (Vect3D), 1, fpio);
    istat += fread(&integral, sizeof (double), 1, fpio);

    if (istat < 7)
        return (NULL);

    tree = newTree3D(data_code, numx, numy, numz, orig.x, orig.y, orig.z, ds.x, ds.y, ds.z, -1.0, integral, NULL);

    istat_cum = 0;
    for (ix = 0; ix < tree->numx; ix++) {
        for (iy = 0; iy < tree->numy; iy++) {
            for (iz = 0; iz < tree->numz; iz++) {
                istat += readNode(fpio, tree->nodeArray[ix][iy][iz]);
                if (istat < 0)
                    return (NULL);
                istat_cum += istat;
            }
        }
    }

    return (tree);

}

/*** function to read an OctNode and all its child nodes ***/

int readNode(FILE *fpio, OctNode* node) {

    int istat;
    int istat_cum;
    int ix, iy, iz;

    float value;


    istat = fread(&value, sizeof (float), 1, fpio); /* node value */
    node->value = (double) value;
    istat += fread(&(node->isLeaf), sizeof (char), 1, fpio); /* leaf flag, 1=leaf */
    if (node->parent != NULL)
        node->level = node->parent->level + 1;
    else
        node->level = 0;

    if (istat < 2)
        return (-1);

    if (node->isLeaf)
        return (1);

    subdivide(node, -1.0, NULL);

    istat_cum = 1;

    for (ix = 0; ix < 2; ix++) {
        for (iy = 0; iy < 2; iy++) {
            for (iz = 0; iz < 2; iz++) {
                if (node->child[ix][iy][iz] != NULL) {
                    istat = readNode(fpio, node->child[ix][iy][iz]);
                    if (istat < 0)
                        return (-1);
                    istat_cum += istat;
                }
            }
        }
    }

    return (istat_cum);

}

/*** function to write a Tree3D ***/

// TODO: modify to support spherical Tree3D!

int writeTree3D(FILE *fpio, Tree3D* tree) {

    int istat;
    int istat_cum;
    int ix, iy, iz;

    istat = fwrite(&(tree->data_code), sizeof (int), 1, fpio);
    istat += fwrite(&(tree->numx), sizeof (int), 1, fpio);
    istat += fwrite(&(tree->numy), sizeof (int), 1, fpio);
    istat += fwrite(&(tree->numz), sizeof (int), 1, fpio);
    istat += fwrite(&(tree->orig), sizeof (Vect3D), 1, fpio);
    istat += fwrite(&(tree->ds), sizeof (Vect3D), 1, fpio);
    istat += fwrite(&(tree->integral), sizeof (double), 1, fpio);

    if (istat < 6)
        return (-1);

    istat_cum = 0;
    for (ix = 0; ix < tree->numx; ix++) {
        for (iy = 0; iy < tree->numy; iy++) {
            for (iz = 0; iz < tree->numz; iz++) {
                istat = writeNode(fpio, tree->nodeArray[ix][iy][iz]);
                if (istat < 0)
                    return (-1);
                istat_cum += istat;
            }
        }
    }

    return (istat_cum);

}

/*** function to write an OctNode and all its child nodes ***/

int writeNode(FILE *fpio, OctNode* node) {

    int istat;
    int istat_cum;
    int ix, iy, iz;

    float value;

    value = (float) node->value;
    istat = fwrite(&value, sizeof (float), 1, fpio); /* node value */
    istat += fwrite(&(node->isLeaf), sizeof (char), 1, fpio); /* leaf flag, 1=leaf */

    if (istat < 2)
        return (-1);

    if (node->isLeaf)
        return (1);

    istat_cum = 1;
    for (ix = 0; ix < 2; ix++) {
        for (iy = 0; iy < 2; iy++) {
            for (iz = 0; iz < 2; iz++) {
                if (node->child[ix][iy][iz] != NULL) {
                    istat = writeNode(fpio, node->child[ix][iy][iz]);
                    if (istat < 0)
                        return (-1);
                    istat_cum += istat;
                }
            }
        }
    }

    return (istat_cum);

}

/** function to get scatter sample for all leafs in results tree
 *  *poct_tree_scatter_volume is weighted volume of cells in results tree
 *  *poct_tree_scatter_volume = SUM(cell_volume * cell_prob / oct_node_value_ref)
 */

int getScatterSampleResultTreeAtLevels(ResultTreeNode* prtree, int value_type, int num_scatter,
        double integral, float* fdata, int npoints, int* pfdata_index,
        double oct_node_value_ref, double *poct_tree_scatter_volume, int level_min, int level_max) {

    OctNode* pnode;
    double xnpoints = 0.0;
    double xval, yval, zval;
    double dx, dy, dz;
    int isample_taken;

    if (prtree->right != NULL)
        npoints = getScatterSampleResultTreeAtLevels(prtree->right, value_type, num_scatter, integral,
            fdata, npoints, pfdata_index, oct_node_value_ref, poct_tree_scatter_volume, level_min, level_max);

    pnode = prtree->pnode;
    //printf("npoints < num_scatter %d && pnode->isLeaf %d && pnode->level >= level_min %d && pnode->level <= level_max %d (level %d min %d max %d\n",
    //        npoints < num_scatter, pnode->isLeaf, pnode->level >= level_min, pnode->level <= level_max, pnode->level, level_min, level_max);
    if (npoints < num_scatter && pnode->isLeaf && pnode->level >= level_min && pnode->level <= level_max) {

        // AJL 20061023 bug fix - prtree value may not be same as pnode value * volume
        //xnpoints = (double) num_scatter * exp(prtree->value - oct_node_value_ref) / integral;
        if (value_type == VALUE_IS_LOG_PROB_DENSITY_IN_NODE) {
            xnpoints = (double) num_scatter * (exp(pnode->value - oct_node_value_ref) * prtree->volume) / integral;
        } else if (value_type == VALUE_IS_PROB_DENSITY_IN_NODE) {
            xnpoints = prtree->volume * (double) num_scatter * (pnode->value / oct_node_value_ref) / integral;
        } else if (value_type == VALUE_IS_PROBABILITY_IN_NODE) {
            // 20140220 AJL - bug fix
            xnpoints = (double) num_scatter * (pnode->value / oct_node_value_ref) / integral;
            //xnpoints = (double) num_scatter * (pnode->value - oct_node_value_ref) / integral;
        }
        //printf("xnpoints %g  num_scatter %d  value %lf  integral %g\n", xnpoints, num_scatter, pnode->value, integral);

        xval = pnode->center.x;
        yval = pnode->center.y;
        zval = pnode->center.z;
        dx = pnode->ds.x / 2.0;
        dy = pnode->ds.y / 2.0;
        dz = pnode->ds.z / 2.0;

        isample_taken = 0;

        //while (xnpoints > 0.0 /*&& npoints < num_scatter*/) {
        while (xnpoints > 0.0 && npoints < num_scatter) { // 20110118 AJL

            if (xnpoints > 1.0 || xnpoints - (double) ((int) xnpoints) > get_rand_double(0.0, 1.0)) {
                fdata[*pfdata_index + 0] = xval + get_rand_double(-dx, dx);
                //printf("npoints %d  *pfdata_index %d  %lf  dx %lf  exp(prtree->value) %le  integral %le\n", npoints, *pfdata_index, fdata[*pfdata_index + 0], dx, exp(prtree->value), integral);
                fdata[*pfdata_index + 1] = yval + get_rand_double(-dy, dy);
                fdata[*pfdata_index + 2] = zval + get_rand_double(-dz, dz);
                fdata[*pfdata_index + 3] = pnode->value;
                //printf("npoints %d  *pfdata_index %d  value %lf  dx %g dy %g dz %g   x %g y %g z %g\n", npoints, *pfdata_index, pnode->value, dx, dy, dz, xval, yval, zval);
                npoints++;
                isample_taken = 1;
                *pfdata_index += 4;
            }

            xnpoints -= 1.0;

        }

        // update weighted scatter volume
        // pnode->value is probability density
        // oct_node_value_ref is maximum probability density
        if (value_type == VALUE_IS_LOG_PROB_DENSITY_IN_NODE) {
            *poct_tree_scatter_volume += prtree->volume * exp(pnode->value - oct_node_value_ref);
        } else if (value_type == VALUE_IS_PROB_DENSITY_IN_NODE) {
            *poct_tree_scatter_volume += prtree->volume * ((pnode->value / oct_node_value_ref) > 0.0 ? (pnode->value / oct_node_value_ref) : 0.0);
        } else if (value_type == VALUE_IS_PROBABILITY_IN_NODE) {
            // 20150215 AJL - bug fix? - need to weight volume by something, TODO: though this (wt by P) is not equivalent to weighting above (wt by PDF)
            *poct_tree_scatter_volume += prtree->volume * ((pnode->value / oct_node_value_ref) > 0.0 ? (pnode->value / oct_node_value_ref) : 0.0);
            // 20140220 AJL - bug fix
            //printf("poct_tree_scatter_volume %f  volume %f  value%f  oct_node_value_ref%f\n", *poct_tree_scatter_volume, prtree->volume, pnode->value, oct_node_value_ref);
            //*poct_tree_scatter_volume += (pnode->value - oct_node_value_ref) > 0.0 ? (pnode->value - oct_node_value_ref) : 0.0;
        }
    }

    if (prtree->left != NULL)
        npoints = getScatterSampleResultTreeAtLevels(prtree->left, value_type, num_scatter, integral,
            fdata, npoints, pfdata_index, oct_node_value_ref, poct_tree_scatter_volume, level_min, level_max);

    return (npoints);
}

/** function to get scatter sample for all leafs in results tree */

int getScatterSampleResultTree(ResultTreeNode* prtree, int value_type, int num_scatter,
        double integral, float* fdata, int npoints, int* pfdata_index,
        double oct_node_value_ref, double *poct_tree_scatter_volume) {

    int level_min = -1;
    int level_max = 9999;

    return (getScatterSampleResultTreeAtLevels(prtree, value_type, num_scatter,
            integral, fdata, npoints, pfdata_index,
            oct_node_value_ref, poct_tree_scatter_volume, level_min, level_max));

}

/** function to integrate exp(val) * volume of all leafs in results tree */

double integrateResultTreeAtLevels(ResultTreeNode* prtree, int value_type, double sum, double oct_node_value_ref, int level_min, int level_max) {

    OctNode* pnode;

    if (prtree->left != NULL)
        sum = integrateResultTreeAtLevels(prtree->left, value_type, sum, oct_node_value_ref, level_min, level_max);

    pnode = prtree->pnode;
    if (pnode->isLeaf && pnode->level >= level_min && pnode->level <= level_max) {
        // result tree value is log(value + volume)
        // AJL 20061023 bug fix - prtree value may not be same as pnode value * volume
        //sum += exp(prtree->value - oct_node_value_ref);
        if (value_type == VALUE_IS_LOG_PROB_DENSITY_IN_NODE) {
            sum += exp(pnode->value - oct_node_value_ref) * prtree->volume;
        } else if (value_type == VALUE_IS_PROB_DENSITY_IN_NODE) {
            sum += prtree->volume * ((pnode->value / oct_node_value_ref) > 0.0 ? (pnode->value / oct_node_value_ref) : 0.0);
        } else if (value_type == VALUE_IS_PROBABILITY_IN_NODE) {
            // 20140220 AJL - bug fix
            sum += (pnode->value / oct_node_value_ref) > 0.0 ? (pnode->value / oct_node_value_ref) : 0.0;
            //sum += (pnode->value - oct_node_value_ref) > 0.0 ? (pnode->value - oct_node_value_ref) : 0.0;
        }
    }

    if (prtree->right != NULL)
        sum = integrateResultTreeAtLevels(prtree->right, value_type, sum, oct_node_value_ref, level_min, level_max);

    return (sum);
}

/** function to integrate exp(val) * volume of all leafs in results tree */

double integrateResultTree(ResultTreeNode* prtree, int value_type, double sum, double oct_node_value_ref) {

    int level_min = -1;
    int level_max = 9999;

    return (integrateResultTreeAtLevels(prtree, value_type, sum, oct_node_value_ref, level_min, level_max));
}

/** function to normalize value of all leafs in results tree to range 0.0-1.0 */

double convertOcttreeValuesToProb(ResultTreeNode* prtree, double integral, double oct_node_value_ref) {

    OctNode* pnode;

    if (prtree->left != NULL)
        integral = convertOcttreeValuesToProb(prtree->left, integral, oct_node_value_ref);

    pnode = prtree->pnode;
    if (pnode->isLeaf) {
        pnode->value = exp(pnode->value - oct_node_value_ref); // replace leaf value with relative prob density
        integral += pnode->value * prtree->volume; // integrate value * cell volume
    }

    if (prtree->right != NULL)
        integral = convertOcttreeValuesToProb(prtree->right, integral, oct_node_value_ref);

    return (integral);

}

/** function to create a new ResultTree using node values */

ResultTreeNode * createResultTree(ResultTreeNode* prtree, ResultTreeNode * pnew_rtree) {

    OctNode* pnode;

    if (prtree->left != NULL)
        pnew_rtree = createResultTree(prtree->left, pnew_rtree);

    pnode = prtree->pnode;
    if (pnode->isLeaf) {
        pnew_rtree = addResult(pnew_rtree, pnode->value, prtree->volume, pnode);
    }

    if (prtree->right != NULL)
        pnew_rtree = createResultTree(prtree->right, pnew_rtree);

    return (pnew_rtree);

}


