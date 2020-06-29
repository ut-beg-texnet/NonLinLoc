// paths

#define PATH1 "/usr/local/lib/fereg96.poly"
#define PATH2 "/usr/local/lib/ch_reg.poly"
#define PATH3 "/usr/local/lib/regions.bin"
#define PATH4 "/usr/local/lib/SED.CAL2"





/* structure for region number and name
   needed by get_region_name_nr  */
typedef struct {
   int start;
   int npairs;
   int reg_number;
   char attribute[12];
   char name[80];
} Aname;

/* globals needed by get_region_name_nr & co  */
Aname *a_name;
int reg_index[1500];
float *la_tab;
float *lo_tab;
int nseg;
int np;

/* function declaration */
int get_region_name_nr(float, float, int, int*, char*, float*, float*);
int inside (float, float, float*, float*, float);
int KSICR(float, float, float, float);
int set_up_names(char*, char*, char*);
Aname region_attributes(float, float);
void celleb(float, float, float*, float*);
void cebell(float, float, float*, float*);
//int csplitstring(char*, char*, char*, char);
