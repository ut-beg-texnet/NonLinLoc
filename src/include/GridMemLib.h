/*------------------------------------------------------------*/
/** grid memory management routines */

typedef struct gridMem {	/* 3D grid data in memory */

	GridDesc* pgrid;	/* pointer to copy of grid description structure */
	void* buffer;		/* corresponding buffer (contiguous floats) */
	void*** array;		/* corresponding array access to buffer */
	int grid_read;		/* grid read flag  = 1 if grid has been read from disk */
	int active;		/* active flag  = 1 if grid is being used in current location */


} GridMemStruct;

/* array of gridMem structure pointers for storing list of grids in memory */
extern GridMemStruct** GridMemList;
extern int GridMemListSize;
extern int GridMemListNumElements;
extern int Num3DGridReadToMemory, MaxNum3DGridMemory;
extern int GridMemListTotalNumElementsAdded;

/* GridLib wrapper functions */
void* NLL_AllocateGrid(GridDesc* pgrid);
void NLL_FreeGrid(GridDesc* pgrid);
void NLL_FreeGridMemory();
void*** NLL_CreateGridArray(GridDesc* pgrid);
void NLL_DestroyGridArray(GridDesc* pgrid);
int NLL_ReadGrid3dBuf(GridDesc* pgrid, FILE* fpio);
GridMemStruct* GridMemList_AddGridDesc(GridDesc* pgrid);
void GridMemList_AddElement(GridMemStruct* pnewGridMemStruct);
void GridMemList_RemoveElementAt(int index);
GridMemStruct* GridMemList_TryToReplaceElementAt(GridMemStruct* pGridMemStruct, GridDesc* pgrid);
GridMemStruct* GridMemList_ElementAt(int index);
int GridMemList_IndexOfGridDesc(int verbose, GridDesc* pgrid);
int GridMemList_NumElements();


/** end of grid memory management routines */
/*------------------------------------------------------------*/


