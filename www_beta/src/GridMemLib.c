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


/*   GridMemLib.c

        Library to manage memory handling of 3D Grid files

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/




#include "GridLib.h"
//#include "ran1.h"
#include "GridMemLib.h"



/*------------------------------------------------------------/ */
/** 3D grid memory management routines to allow persistence of grids in memory */

#define USE_GRID_LIST 1
#define GRIDMEM_MESSAGE 2

/*** wrapper function to allocate buffer for 3D grid ***/

void* NLL_AllocateGrid(GridDesc* pgrid) {
    int index, nactive, ngrid_read, n;
    void* fptr = NULL;
    GridMemStruct* pGridMemStruct = NULL;

    //printf("IN: NLL_AllocateGrid\n");

    if (USE_GRID_LIST) {

        if ((index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
            // already in list
            pGridMemStruct = GridMemList_ElementAt(index);
            pGridMemStruct->active = 1;
            fptr = pGridMemStruct->buffer;
            if (message_flag >= GRIDMEM_MESSAGE)
                printf("GridMemManager: Grid exists in mem (%d/%d): %s\n", index, GridMemListNumElements, pGridMemStruct->pgrid->title);
            return (fptr);
        } else {
            // check number of active grids in list
            nactive = 0;
            ngrid_read = 0;
            for (n = 0; n < GridMemList_NumElements(); n++) {
                pGridMemStruct = GridMemList_ElementAt(n);
                nactive += pGridMemStruct->active;
                ngrid_read += pGridMemStruct->grid_read;
            }
            // list already full of active grids, do normal allocation
            if (MaxNum3DGridMemory > 0 && nactive >= MaxNum3DGridMemory) {
                fptr = AllocateGrid(pgrid);
                if (message_flag >= GRIDMEM_MESSAGE)
                    printf("GridMemManager: Memory full (%d/%d): %s\n", index, GridMemListNumElements, pGridMemStruct->pgrid->title);
                return (fptr);
            }
            // ok
            // remove an inactive grid if necessary
            if (MaxNum3DGridMemory > 0 && ngrid_read >= MaxNum3DGridMemory) {
                for (n = GridMemList_NumElements() - 1; n >= 0; n--) {
                    pGridMemStruct = GridMemList_ElementAt(n);
                    if (!pGridMemStruct->active && pGridMemStruct->grid_read) {
                        GridMemList_RemoveElementAt(n);
                        break;
                    }
                }
            }
            // create new list element
            pGridMemStruct = GridMemList_AddGridDesc(pgrid);
            fptr = pGridMemStruct->buffer;
            if (fptr == NULL) {
                // error allocating grid memory or out of memory
                GridMemList_RemoveElementAt(GridMemList_NumElements() - 1);
            }
            return (fptr);
        }

    } else {
        fptr = AllocateGrid(pgrid);
        return (fptr);
    }

}

/*** wrapper function to free buffer for 3D grid ***/

void NLL_FreeGrid(GridDesc* pgrid) {
    int index;
    GridMemStruct* pGridMemStruct;

    //printf("IN: NLL_FreeGrid\n");
    if (USE_GRID_LIST && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
        pGridMemStruct = GridMemList_ElementAt(index);
        pGridMemStruct->active = 0;
        //pgrid->buffer = NULL;
        return;
    }

    FreeGrid(pgrid);
}


/*** return element from GridMemList ***/
// 20131227 AJL - bug fix, added this function.  Before, grid memory was not freed.

void NLL_FreeGridMemory() {

    int index;
    int numElements = GridMemListNumElements;
    // 20141219 AJL - bug fix, GridMemListNumElements is decremented each time GridMemList_RemoveElementAt() called!)
    //for (index = 0; index < GridMemListNumElements; index++) {
    for (index = numElements; index >= 0; index--) {  // 20141219 AJL - bug fix, must remove from end
        GridMemList_RemoveElementAt(index);
    }
    free(GridMemList);  // 20141219 AJL - bug fix, added this free

}

/*** wrapper function to create array for accessing 3D grid ***/

void*** NLL_CreateGridArray(GridDesc* pgrid) {

    void*** fptr = NULL;

    int index;
    GridMemStruct* pGridMemStruct;

    //printf("IN: NLL_CreateGridArray\n");
    if (USE_GRID_LIST && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
        pGridMemStruct = GridMemList_ElementAt(index);
        fptr = pGridMemStruct->array;
    } else
        fptr = CreateGridArray(pgrid);

    return (fptr);
}

/*** wrapper function to free array for accessing 3D grid ***/

void NLL_DestroyGridArray(GridDesc* pgrid) {

    int index;

    //printf("IN: NLL_DestroyGridArray\n");
    if (USE_GRID_LIST && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
        pgrid->array = NULL;
        return;
    }

    DestroyGridArray(pgrid);
}

/*** wrapper function to read entire grid buffer from disk ***/

int NLL_ReadGrid3dBuf(GridDesc* pgrid, FILE* fpio) {

    int istat;

    int index;
    GridMemStruct* pGridMemStruct;

    //printf("IN: NLL_ReadGrid3dBuf\n");
    if (USE_GRID_LIST && (index = GridMemList_IndexOfGridDesc(0, pgrid)) >= 0) {
        pGridMemStruct = GridMemList_ElementAt(index);
        if (!pGridMemStruct->grid_read) {
            istat = ReadGrid3dBuf(pGridMemStruct->pgrid, fpio);
            pGridMemStruct->grid_read = 1;
        }
    } else {
        istat = ReadGrid3dBuf(pgrid, fpio);
    }

    return (0);
}

/*** add GridDescription to GridMemList ***/

GridMemStruct* GridMemList_AddGridDesc(GridDesc* pgrid) {
    GridMemStruct* pnewGridMemStruct;

    //printf("IN: GridMemList_AddGridDesc\n");
    pnewGridMemStruct = (GridMemStruct*) malloc(sizeof (GridMemStruct));
    pnewGridMemStruct->pgrid = (GridDesc*) malloc(sizeof (GridDesc));
    *(pnewGridMemStruct->pgrid) = *pgrid;
    strcpy(pnewGridMemStruct->pgrid->chr_type, pgrid->chr_type);
    strcpy(pnewGridMemStruct->pgrid->title, pgrid->title);
    pnewGridMemStruct->buffer = AllocateGrid(pnewGridMemStruct->pgrid);
    pnewGridMemStruct->array = CreateGridArray(pnewGridMemStruct->pgrid);
    pnewGridMemStruct->active = 1;
    pnewGridMemStruct->grid_read = 0;

    GridMemList_AddElement(pnewGridMemStruct);

    return (pnewGridMemStruct);

}



/*** add element to GridMemList ***/

#define LIST_SIZE_INCREMENT 10

void GridMemList_AddElement(GridMemStruct* pnewGridMemStruct) {
    int n;
    int newGridMemListSize;
    GridMemStruct** newGridMemList;

    //printf("IN: GridMemList_AddElement\n");
    if (GridMemListSize <= GridMemListNumElements) {
        // allocate enlarged list
        newGridMemListSize = GridMemListSize + LIST_SIZE_INCREMENT;
        newGridMemList = (GridMemStruct**)
                malloc(newGridMemListSize * sizeof (GridMemStruct*));
        // load old list to new list
        for (n = 0; n < GridMemListSize; n++)
            newGridMemList[n] = GridMemList[n];
        for (n = GridMemListSize; n < newGridMemListSize; n++)
            newGridMemList[n] = NULL;
        GridMemListSize = newGridMemListSize;
        if (GridMemList != NULL)
            free(GridMemList);
        GridMemList = newGridMemList;
    }

    // load new element
    GridMemList[GridMemListNumElements] = pnewGridMemStruct;
    GridMemListNumElements++;

    if (message_flag >= GRIDMEM_MESSAGE)
        printf("GridMemManager: Add grid (%d): %s\n", GridMemListNumElements - 1, pnewGridMemStruct->pgrid->title);
}

/*** remove element from GridMemList ***/

void GridMemList_RemoveElementAt(int index) {
    int n;
    GridMemStruct* pGridMemStruct;

    //printf("IN: GridMemList_RemoveElementAt\n");
    // check ranges
    if (index < 0 || index >= GridMemListNumElements)
        return;

    // free allocated memory
    pGridMemStruct = GridMemList[index];
    if (message_flag >= GRIDMEM_MESSAGE)
        printf("GridMemManager: Remove grid (%d/%d): %s\n", index, GridMemListNumElements, pGridMemStruct->pgrid->title);
    DestroyGridArray(pGridMemStruct->pgrid);
    FreeGrid(pGridMemStruct->pgrid);
    free(pGridMemStruct->pgrid);
    free(pGridMemStruct);


    // shift down element references
    for (n = index; n < GridMemListNumElements - 1; n++)
        GridMemList[n] = GridMemList[n + 1];

    GridMemList[n] = NULL;
    GridMemListNumElements--;

    return;
}

/*** return element from GridMemList ***/

GridMemStruct* GridMemList_ElementAt(int index) {

    //printf("IN: GridMemList_ElementAt\n");
    // check ranges
    if (index < 0 || index >= GridMemListNumElements)
        return (NULL);

    return (GridMemList[index]);
}

/*** find index of grid desc in GridMemList ***/

int GridMemList_IndexOfGridDesc(int verbose, GridDesc* pgrid) {

    int n;

    //printf("IN: GridMemList_IndexOfGridDesc\n");
    for (n = 0; n < GridMemListNumElements; n++) {
        if (verbose) printf("indexOf: %s ==? %s\n", GridMemList[n]->pgrid->title, pgrid->title);
        if (strcmp(GridMemList[n]->pgrid->title, pgrid->title) == 0)
            return (n);
    }

    if (verbose) printf("indexOf: NOT FOUND\n");
    return (-1);

}

/*** return size of GridMemList ***/

int GridMemList_NumElements() {

    //printf("IN: GridMemList_NumElements\n");
    return (GridMemListNumElements);

}



/** end of 3D grid memory management routines */
/*------------------------------------------------------------/ */


