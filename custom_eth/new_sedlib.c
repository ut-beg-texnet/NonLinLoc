/*================================================
  Library for reading KP- and GSE-files
 ================================================*/
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <string.h>

#ifdef Motif
#define free XtFree
#define malloc XtMalloc
#define realloc XtRealloc
#endif

#include "new_sedlib.h"

#define DEBUG 0

/*=================================================================*/

int NewOpenFile(filename)
char *filename;
{

    int i;
    unsigned char *file_data;
    int file_length;
    char filname[256];

    PointerInitialize();

    strcpy(filname, sedlibPath);
    if (strlen(filname) < 1) {
        strcpy(filname, "/app/lib/");
    }
    strcat(filname, "extensions.def");
    nfext = ReadExtensions(filname);

    strcpy(workingfile, filename);

    /* check for compressed file: */
    CeckCompressed(workingfile, &is_compressed);

    /* check for extension file: */
    LastExt(workingfile, cleanName, extension);

    if (!((i = strlen(extension)) == 0 || (strncmp(extension, "GSE", 3) == 0))) {
        if (is_compressed) {
            unlink(workingfile);
        }
        return (1);
    }
    if (strlen(extension) == 0) {
        /* read full data file into memory */
        UpdatePermission = GetPermission(workingfile, "file");
        i = ReadFullFile(workingfile, &file_data, &file_length);
        if (i > 0) {
            i = fprintf(stderr, "Error: ReadFullFile returned with %d\n", i);
            return (2);
        }
        i = DecodeKPfile(file_data, file_length);
        if (i > 0) {
            i = fprintf(stderr, "Error: DecodeKPfile returned with %d\n", i);
            return (2);
        }
        free(file_data);
        IsGSEFormat = False;
    } else {
        UpdatePermission = GetPermission(workingfile, "dir");
        i = DecodeGSEfile(workingfile, cleanName);
        if (i > 0) {
            i = fprintf(stderr, "Error: DecodeGSEfile returned with %d\n", i);
            return (2);
        }
        IsGSEFormat = True;
    }
    i = AddRecSys("instrument_type.lib");
    strcpy(filname, sedlibPath);
    strcat(filname, "SED.CAL2");
    i = ReadCALfromFile(filname, "", 0);
    i = AddCalibToSignals();
    i = AddRecSysToPicks();
    return (0);
}

BoolSED GetPermission(file, what)
char *file;
char *what;
{
    char fcpy[256], dn[256];
    int i = 0;
    int mode;
    BoolSED Permission;

    if (strcmp(what, "file") == 0) {
        mode = R_OK | W_OK;
        i = access(file, R_OK | W_OK);
    }
    if (strcmp(what, "dir") == 0) {
        strcpy(fcpy, file);
        strcpy(dn, dirname(fcpy));
        mode = R_OK | W_OK | X_OK;
        i = access(dn, mode);
    }

    if (i == 0) Permission = True;
    else Permission = False;
    return ( Permission);
}

int ReadExtensions(file)
char *file;
{
    FILE *Fop;
    int limit, n, i;

    if (nfext > 0) return (nfext);
    if ((Fop = fopen(file, "r")) == NULL) { /* set to defaults */
        n = fprintf(stderr, "file %s not found, default extensions used instead\n", file);
        strcpy(fext[TRIGGER].extension, "TRIGGER");
        strcpy(fext[APICK].extension, "AUTOPICK");
        strcpy(fext[MPICK].extension, "MANUPICK");
        strcpy(fext[AFILTER].extension, "AUTOFILTER");
        strcpy(fext[MFILTER].extension, "MANUFILTER");
        strcpy(fext[ALOC].extension, "AUTOLOC");
        strcpy(fext[APDE].extension, "AUTOPDE");
        strcpy(fext[MLOC].extension, "MANULOC");
        strcpy(fext[MPDE].extension, "MANUPDE");
        strcpy(fext[MAGNITUDE].extension, "MAGNITUDE");
        strcpy(fext[TIMESTATUS].extension, "TIMESTATUS");
        strcpy(fext[SIGNALSTATUS].extension, "SIGNALSTATUS");
        strcpy(fext[WAVEFORM].extension, "GSE2");
        strcpy(fext[NOEXT].extension, "");
        fext[TRIGGER].limit = 0;
        fext[APICK].limit = 0;
        fext[MPICK].limit = 0;
        fext[AFILTER].limit = 0;
        fext[MFILTER].limit = 0;
        fext[ALOC].limit = 0;
        fext[APDE].limit = 0;
        fext[MLOC].limit = 0;
        fext[MPDE].limit = 0;
        fext[MAGNITUDE].limit = 0;
        fext[TIMESTATUS].limit = 0;
        fext[SIGNALSTATUS].limit = 0;
        fext[NOEXT].limit = 0;
        return (NOEXT);
    }
    fgets(line, 256, Fop);
    while (feof(Fop) == 0) {
        for (i = 0; i < strlen(line); i++) {
            if (line[i] == '\t') line[i] = ' ';
        }
        n = sscanf(line, "%s %s %d", str1, str2, &limit);
        if (strcmp(str1, "triggers") == 0) {
            strcpy(fext[TRIGGER].extension, str2);
            fext[TRIGGER].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "autopicks") == 0) {
            strcpy(fext[APICK].extension, str2);
            fext[APICK].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "manualpicks") == 0) {
            strcpy(fext[MPICK].extension, str2);
            fext[MPICK].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "autofilter") == 0) {
            strcpy(fext[AFILTER].extension, str2);
            fext[AFILTER].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "manualfilter") == 0) {
            strcpy(fext[MFILTER].extension, str2);
            fext[MFILTER].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "autolocation") == 0) {
            strcpy(fext[ALOC].extension, str2);
            fext[ALOC].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "autoloclist") == 0) {
            strcpy(fext[APDE].extension, str2);
            fext[APDE].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "manuallocation") == 0) {
            strcpy(fext[MLOC].extension, str2);
            fext[MLOC].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "manualloclist") == 0) {
            strcpy(fext[MPDE].extension, str2);
            fext[MPDE].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "magnitude") == 0) {
            strcpy(fext[MAGNITUDE].extension, str2);
            fext[MAGNITUDE].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "timestatus") == 0) {
            strcpy(fext[TIMESTATUS].extension, str2);
            fext[TIMESTATUS].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "signalstatus") == 0) {
            strcpy(fext[SIGNALSTATUS].extension, str2);
            fext[SIGNALSTATUS].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
        if (strcmp(str1, "waveform") == 0) {
            strcpy(fext[WAVEFORM].extension, str2);
            fext[WAVEFORM].limit = limit;
            fgets(line, 256, Fop);
            continue;
        }
    }
    strcpy(fext[NOEXT].extension, "");
    fext[NOEXT].limit = 0;
    fclose(Fop);
    return (NOEXT);
}

int AddComp2Table(comp)
char *comp;
{
    int comp_bit, k;

    /* component table */
    comp[3] = 0;
    /* put this component into table if it is not yet there */
    comp_bit = -1;
    for (k = 0; k < ncomp_table; k++) {
        if (strcmp(comp, comp_table[k]) == 0) {
            comp_bit = k;
            break;
        }
    }
    if (comp_bit < 0) {
        if (ncomp_table < 32) { /* not yet in table */
            strcpy(comp_table[ncomp_table], comp);
            comp_bit = ncomp_table;
            ncomp_table++;
        } else comp_bit = 31;
    }

    return ( comp_bit);
}

int ReadFullFile(workingfile, file_d, file_length)
char *workingfile;
unsigned char **file_d;
int *file_length;
{
    FILE *fp_data = NULL;
    unsigned char *file_data;
    struct stat statbuf;
    int file_l;
    int nbytes, i;
    *file_length = 0;

    /* get file length */
    if (stat(workingfile, &statbuf) == 0) {
        file_l = statbuf.st_size;
    } else {
        return (1);
    }
    /*
     * check for ordinary file type:
     */
    if (!S_ISREG(statbuf.st_mode)) {
        return (2);
    }
    /**
     ** now open and read the data file
     **/
    if ((fp_data = fopen(workingfile, "r")) == NULL) {
        return (3);
    }

    /*
     * allocat memory for the whole file:
     */
    file_data = (unsigned char *) malloc(file_l);
    if (file_data == NULL) {
        i = fprintf(stderr, "Error: not enough memory for data.\n");
        return (4);
    };
    /*
     * read the data file:
     */
    nbytes = fread(file_data, sizeof (char), file_l, fp_data);
    if (nbytes != file_l) {
        i = fprintf(stderr, "Error: file length is %d, bytes read are %d\n", file_l, nbytes);
        return (5);
    }
    *file_length = nbytes;
    *file_d = file_data;
    return (0);
}

void CeckCompressed(filename, is_comp)
char *filename;
BoolSED *is_comp;
{
    char dirn[256], basen[256], fclean[256], extension[256];
    char decompress[256], runstring[256], fn[256];
    int i;
    BoolSED compressed;

    strcpy(fn, filename);
    strcpy(basen, basename(fn));
    strcpy(fn, filename);
    strcpy(dirn, dirname(fn));
    LastExt(basen, fclean, extension);
    /** treat compressed data files: extension .Z .gz .zip **/
    decompress[0] = 0;
    if (strcmp(extension, "Z") == 0) {
        strcpy(decompress, "uncompress");
    }
    if (strcmp(extension, "gz") == 0) {
        strcpy(decompress, "gunzip");
    }
    if (strcmp(extension, "zip") == 0) {
        strcpy(decompress, "unzip -d /tmp");
    }
    if ((i = strlen(decompress)) > 0) {
        strcpy(dirn, "/tmp/");
        strcat(dirn, basen);
        i = sprintf(runstring, "cp %s %s;%s %s",
                filename, dirn, decompress, dirn);
        i = system(runstring);
        strcpy(dirn, "/tmp/");
        strcat(dirn, fclean);
        strcpy(filename, dirn);
        compressed = True;
    } else {
        compressed = False;
    }
    *is_comp = compressed;
}

void LastExt(string, clean, extens)
char *string, *clean, *extens;
{
    int i, j;
    char dname[256], bname[256], fname[256];

    strcpy(fname, string);
    strcpy(dname, dirname(fname));
    strcpy(fname, string);
    strcpy(bname, basename(fname));

    j = strlen(bname) - 1;
    strcpy(fname, bname);
    extens[0] = 0;

    for (i = j; i > 0; i--) {
        if (bname[i] == '.') {
            strcpy(extens, &bname[i + 1]);
            fname[i] = 0;
            break;
        }
    }
    for (i = 0; i < strlen(extens); i++) {
        j = toupper(extens[i]);
        extens[i] = (char) j;
    }
    strcpy(clean, dname);
    strcat(clean, "/");
    strcat(clean, fname);
}

void PointerInitialize() {
    int i, j;
    if (nmxTrigger != NULL) free(nmxTrigger);
    nmxTrigger = NULL;
    NnmxTrigger = 0;
    if (autoPicks != NULL) free(autoPicks);
    autoPicks = NULL;
    NautoPicks = 0;
    if (manualPicks != NULL) free(manualPicks);
    manualPicks = NULL;
    NmanualPicks = 0;
    if (autofil != NULL) free(autofil);
    autofil = NULL;
    Nautofil = 0;
    if (manufil != NULL) free(manufil);
    manufil = NULL;
    Nmanufil = 0;
    if (calib != NULL) free(calib);
    calib = NULL;
    ncalib = 0;
    for (i = 0; i < nstat; i++) {
        for (j = 0; j < station[i].nsig; j++) {
            if (station[i].wave[j] != NULL) {
                free(station[i].wave[j]);
                station[i].wave[j] = NULL;
            }
            if (station[i].shead[j].status_out != NULL) {
                free(station[i].shead[j].status_out);
                station[i].shead[j].status_out = NULL;
                free(station[i].shead[j].status_back);
                station[i].shead[j].status_back = NULL;
            }
            if (station[i].shead[j].timing_out != NULL) {
                free(station[i].shead[j].timing_out);
                station[i].shead[j].timing_out = NULL;
                free(station[i].shead[j].timing_back);
                station[i].shead[j].timing_back = NULL;
            }
        }
    }
    memset(station, 0, MAX_STAT * sizeof (Stat_Str));
    memset(&locNmxT, 0, sizeof (NewLocPar));
    memset(&locAutoP, 0, sizeof (NewLocPar));
    memset(&locManualP, 0, sizeof (NewLocPar));
    memset(&autoLoc, 0, sizeof (NewLoc));
    memset(&manuLoc, 0, sizeof (NewLoc));
    memset(&mag_fixed, 0, sizeof (mag_fixed));
}

int AddCalibToSignals() {
    int i, j, k, l;
    int overlap, value;
    double t1, t2, c1, c2;

    t1 = (double) time_minimum;
    t2 = (double) time_maximum;
    for (i = 0; i < nstat; i++) {
        for (j = 0; j < station[i].nsig; j++) {
            for (k = 0; k < ncalib; k++) {
                if (strcmp(calib[k].name, station[i].name) != 0) continue;
                if (strcmp(calib[k].comp, station[i].shead[j].s_comp) != 0) continue;
                c1 = (double) calib[k].time1;
                c2 = (double) calib[k].time2;
                overlap = TimeWindowOverLap(t1, t2, c1, c2, &value);
                if (overlap) {
                    station[i].shead[j].calib = &calib[k];
                    for (l = 0; l < NrecSys; l++) {
                        if (calib[k].auxid[0] == recSys[l].recsys[0]) {
                            station[i].shead[j].s_resol = recSys[l].resol;
                            station[i].shead[j].s_saturation = recSys[l].clip;
                            if (station[i].shead[j].s_volts < 0) {
                                station[i].shead[j].s_volts = -recSys[l].volts;
                            } else {
                                station[i].shead[j].s_volts = recSys[l].volts;
                            }
                            strcpy(station[i].shead[j].s_recsys, recSys[l].recsys);
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    return (0);
}

void GetCal4StnTime(stn, time, idx, Nidx)
char *stn;
int time;
int *idx, *Nidx;
{
    int i, k;

    i = 0;
    if (DEBUG) printf("TP GC4ST 1 ncalib=%d\n", ncalib);
    for (k = 0; k < ncalib; k++) {
        //if (DEBUG) printf("TP GC4ST 2 %s %s\n", calib[k].name,stn);
        if (strcmp(calib[k].name, stn) != 0) continue;
        if (DEBUG) printf("TP GC4ST 3 %d %d %d\n", time, calib[k].time1, calib[k].time2);
        if (time < calib[k].time1 || time > calib[k].time2) continue;
        idx[i] = k;
        i++;
    }
    *Nidx = i;
    if (DEBUG) printf("TP GC4ST 4 %d\n", *Nidx);
}

int AddRecSysToPicks() {
    int i;
    for (i = 0; i < NnmxTrigger; i++) {
        strcpy(nmxTrigger[i].recsys, "NAQ");
    }
    i = AddRecSysToThesePicks(autoPicks, NautoPicks);
    i = AddRecSysToThesePicks(manualPicks, NmanualPicks);
    /************************************************
    for (i=0;i<NautoPicks;i++) {
       for (j=0;j<nstat;j++) {
          if( strcmp(station[j].name,autoPicks[i].stname) == 0 ) {
             if ( station[j].shead[0].calib != NULL )
                  strcpy(autoPicks[i].recsys,station[j].shead[0].calib->recsys);
             break;
          }
       }
       if (autoPicks[i].recsys[0] == 0 ) strcpy(autoPicks[i].recsys," ");
    }
    for (i=0;i<NmanualPicks;i++) {
       for (j=0;j<nstat;j++) {
          if( strcmp(station[j].name,manualPicks[i].stname) == 0 ) {
             if ( station[j].shead[0].calib != NULL )
                strcpy(manualPicks[i].recsys,station[j].shead[0].calib->recsys);
             break;
          }
       }
       if (manualPicks[i].recsys[0] == 0 ) strcpy(manualPicks[i].recsys," ");
    }
     ************************************************/
    return (0);
}

int AddRecSysToThesePicks(p, n)
NewPick *p;
int n;
{
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < nstat; j++) {
            if (strcmp(station[j].name, p[i].stname) == 0) {
                if (station[j].shead[0].calib != NULL)
                    strcpy(p[i].recsys, station[j].shead[0].calib->recsys);
                break;
            }
        }
        if (p[i].recsys[0] == 0) strcpy(p[i].recsys, " ");
    }
    return (0);

}

int AddRecSys(file)
char *file;
{
    FILE *fp = NULL;
    int i, l;
    char lfile[256];

    strcpy(lfile, sedlibPath);
    if (strlen(file) == 0) {
        strcat(lfile, "instrument_type.lib");
    } else {
        strcat(lfile, file);
    }
    if (recSys != 0) {
        free(recSys);
        recSys = NULL;
    }
    NrecSys = 0;

    if ((fp = fopen(file, "r")) != NULL) {
        fgets(line, 256, fp);
        while (!feof(fp)) {
            if (line[0] != '#') {
                recSys = (RecSysDefinition *) realloc(recSys, (NrecSys + 1) * sizeof (RecSysDefinition));
                memset(&recSys[NrecSys], 0, sizeof (RecSysDefinition));
                for (i = 0; i < strlen(line); i++) {
                    if (line[i] == '\t') line[i] = ' ';
                }
                l = sscanf(line, "%d %s %f %f", &recSys[NrecSys].resol,
                        recSys[NrecSys].recsys,
                        &recSys[NrecSys].volts,
                        &recSys[NrecSys].clip);
                NrecSys++;
            }
            fgets(line, 256, fp);
        }
        fclose(fp);
    } else {
        NrecSys = 11;
        recSys = (RecSysDefinition *) realloc(recSys, NrecSys * sizeof (RecSysDefinition));
        memset(&recSys[0], 0, NrecSys * sizeof (RecSysDefinition));
        recSys[0].resol = 12;
        recSys[1].resol = 12;
        recSys[2].resol = 12;
        recSys[3].resol = 10;
        recSys[4].resol = 19;
        recSys[5].resol = 24;
        recSys[6].resol = 24;
        recSys[7].resol = 24;
        recSys[8].resol = 24;
        recSys[9].resol = 16;
        recSys[10].resol = 17;
        recSys[0].volts = 2;
        recSys[1].volts = 0.5;
        recSys[2].volts = 10;
        recSys[3].volts = 2;
        recSys[4].volts = 2.5;
        recSys[5].volts = 20;
        recSys[6].volts = 0.5;
        recSys[7].volts = 0.6;
        recSys[8].volts = 2;
        recSys[9].volts = 0.032768;
        recSys[10].volts = 320;
        recSys[0].clip = 0.7;
        recSys[1].clip = 1.0;
        recSys[2].clip = 1.0;
        recSys[3].clip = 0.7;
        recSys[4].clip = 1.0;
        recSys[5].clip = 1.0;
        recSys[6].clip = 1.0;
        recSys[7].clip = 1.0;
        recSys[8].clip = 1.0;
        recSys[9].clip = 1.0;
        recSys[10].clip = 1.0;
        strcpy(recSys[0].recsys, "AHP");
        strcpy(recSys[1].recsys, "MLR");
        strcpy(recSys[2].recsys, "LEE");
        strcpy(recSys[3].recsys, "FHP");
        strcpy(recSys[4].recsys, "PDR");
        strcpy(recSys[5].recsys, "NAQ");
        strcpy(recSys[6].recsys, "NAQ");
        strcpy(recSys[7].recsys, "NAQ");
        strcpy(recSys[8].recsys, "NAQ");
        strcpy(recSys[9].recsys, "M88");
        strcpy(recSys[10].recsys, "LEE7");
    }
    return (0);
}

int TimeWindowOverLap(z1, z2, x1, x2, type)
double z1, z2, x1, x2;
int *type;
{
    /*****************************************
      |=======|        |-------|
      z1      z2       x1      x2
     
      |--------|  |=========|  :
          z1-x1 <0; z1-x2 <0; z2-x1 <0; z2-x2 <0
            -2    +  -1     +   -1    +   -1     = -5

      |========|  |---------|  :
          z1-x1 >0; z1-x2 >0; z2-x1 >0; z2-x2 >0
            +2    +   +1    +   +1    +   +1     = 5
   
      |========|
           |---------|:
          z1-x1 <0; z1-x2 <0; z2-x1 >0; z2-x2 <0
            -2    +   -1    +   +1    +   -1     = -3

      |--------|
           |========|:
          z1-x1 >0; z1-x2 <0; z2-x1 >0; z2-x2 >0
            +2    +   -1    +   +1    +   +1     = 3
     
      |------------|
        |========|:
          z1-x1 >0; z1-x2 <0; z2-x1 >0; z2-x2 <0
            +2    +   -1    +   +1    +   -1     = 1
     
      |============|
        |--------|:
          z1-x1 <0; z1-x2 <0; z2-x1 >0; z2-x2 >0
            -2    +   -1    +   +1    +   +1     = -1
     
     *****************************************/
    int k = 0, ov = 0;
    if (z1 < x1) k -= 2;
    else k += 2;
    if (z1 < x2) k--;
    else k++;
    if (z2 < x1) k--;
    else k++;
    if (z2 < x2) k--;
    else k++;
    if (abs(k) < 5) ov = 1;
    *type = k;
    return (ov);
}

int DefaultLocPar(lpar)
NewLocPar *lpar;
{
    lpar->n_picks = 0;
    lpar->st_list = 1;
    lpar->model = 1;
    lpar->lat = 0;
    lpar->lon = 0;
    lpar->depth = 10;
    lpar->xnear = 200;
    lpar->xfar = 300;
    lpar->pos = 1.71;
    lpar->dray = 0;
    lpar->daz = 0;
    strcpy(lpar->type, "Unkn");
    lpar->fix = 0;
    lpar->yr = 0;
    lpar->mo = 0;
    lpar->dy = 0;
    lpar->hr = 0;
    lpar->mn = 0;
    lpar->sec = 0;
    lpar->swt = 1;
    lpar->AmpFilt = 0;
    lpar->NoMag = 0;
    lpar->dumy2 = 0;

    return (0);
}

/************************************
Decode GSE Files
 *************************************/
int DecodeGSEfile(file, wfile)
char *file, *wfile;
{
    int i;

    if ((i = ReadWIDfromFile(file, "", 0)) > 0) return (i);
    i = ReadCALfromFile(file, "", 0);
    i = ReadPicksfromFile(wfile, fext[TRIGGER].extension, fext[TRIGGER].limit, &nmxTrigger, &NnmxTrigger, &locNmxT);
    i = ReadPicksfromFile(wfile, fext[APICK].extension, fext[APICK].limit, &autoPicks, &NautoPicks, &locAutoP);
    i = ReadPicksfromFile(wfile, fext[MPICK].extension, fext[MPICK].limit, &manualPicks, &NmanualPicks, &locManualP);
    i = ReadFiltersfromFile(wfile, fext[AFILTER].extension, fext[AFILTER].limit, &autofil, &Nautofil);
    i = ReadFiltersfromFile(wfile, fext[MFILTER].extension, fext[MFILTER].limit, &manufil, &Nmanufil);
    i = ReadLocationfromFile(wfile, fext[ALOC].extension, fext[ALOC].limit, &autoLoc);
    i = ReadLocationfromFile(wfile, fext[MLOC].extension, fext[MLOC].limit, &manuLoc);
    i = ReadMagnitudefromFile(wfile, fext[MAGNITUDE].extension, fext[MAGNITUDE].limit, &mag_fixed);
    i = ReadTimeStatusFile(wfile, fext[TIMESTATUS].extension, fext[TIMESTATUS].limit);
    i = ReadSignalStatusFile(wfile, fext[SIGNALSTATUS].extension, fext[SIGNALSTATUS].limit);
    return (0);
}

int ReadMagnitudefromFile(file, typ, limit, m_fix)
char *file, *typ;
int limit;
FixMag *m_fix;
{
    FILE *Fop;
    char mtyp[5];
    int i, k;
    float mag;

    strcpy(str4, file);
    Fop = OpenLimiter(str4, typ, "r", limit);
    if (Fop != NULL) {
        fgets(str2, 256, Fop);
        while (!feof(Fop)) {
            if (str2[0] != '#') {
                str2[54] = 0;
                while (strlen(str2) > 0) {
                    csplitstring(str2, str1, str2, '=');
                    i = sscanf(str2, "%f", &mag);
                    if (i == 1) {
                        k = 0;
                        for (i = 0; i < strlen(str1); i++) {
                            if ((str1[i] > 64 && str1[i] < 91) || (str1[i] > 96 && str1[i] < 123)) {
                                mtyp[k] = str1[i];
                                k++;
                            }
                        }
                        mtyp[k] = 0;
                    } else memset(&mtyp, 0, 5);
                    if (mtyp[1] == 'b' || mtyp[1] == 'B') {
                        if (mag < 10) mag_fixed.fix_mb = mag;
                    } else if (mtyp[1] == 'l' || mtyp[1] == 'L') {
                        if (mag < 10) mag_fixed.fix_ml = mag;
                    } else if (mtyp[0] != 0) {
                        if (mag < 10) {
                            mag_fixed.fix_ms = mag;
                            strcpy(mag_fixed.fix_ms_str, mtyp);
                        }
                    }
                }
            }
            fgets(str2, 256, Fop);
        }
        fclose(Fop);
        return ( 0);
    }
    return ( 1);
}

int ReadTimeStatusFile(file, typ, limit)
char *file, *typ;
int limit;
{
    FILE *Fop;

    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    if (Fop != NULL) {

        fclose(Fop);
        return ( 0);
    }
    return ( 1);
}

int ReadSignalStatusFile(file, typ, limit)
char *file, *typ;
int limit;
{
    FILE *Fop;
    char stn[8], chan[8];
    int i, j, n1, n2;

    strcpy(str2, file);
    Fop = OpenLimiter(str2, typ, "r", limit);
    if (Fop != NULL) {
        fgets(line, 256, Fop);
        while (!feof(Fop)) {
            if (line[0] == '#') {
                fgets(line, 256, Fop);
                continue;
            }
            sscanf(line, "%s %s %d %d", stn, chan, &n1, &n2);
            for (i = 0; i < nstat; i++) {
                if (strcmp(station[i].name, stn) == 0) {
                    for (j = 0; j < station[i].nsig; j++) {
                        if (strcmp(station[i].shead[j].s_comp, chan) == 0) {
                            station[i].shead[j].status_out = (int *) realloc(station[i].shead[j].status_out, sizeof (int) *(station[i].shead[j].Nstatus_out + 1));
                            station[i].shead[j].status_back = (int *) realloc(station[i].shead[j].status_back, sizeof (int) *(station[i].shead[j].Nstatus_out + 1));
                            station[i].shead[j].status_out[station[i].shead[j].Nstatus_out] = n1;
                            station[i].shead[j].status_back[station[i].shead[j].Nstatus_out] = n2;
                            station[i].shead[j].Nstatus_out++;
                            break;
                        }
                    }
                    break;
                }
            }
            fgets(line, 256, Fop);
        }
        fclose(Fop);
        return ( 0);
    }
    return ( 1);
}

int ReadPicksfromFile(file, typ, limit, p, Np, locP)
char *file, *typ;
int limit;
NewPick **p;
int *Np;
NewLocPar *locP;
{
    FILE *Fop;
    char key[5];
    int i, j, l, n, npick, yr, mo, dy, hr, mn, evnr, eminute = 0;
    float tmin;

    NewPick *lp = NULL;

    tmin = 0.2;
    evnr = 0;
    npick = 0;
    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    if (Fop != NULL) {
        locP->AmpFilt = 'T';
        fgets(line, 256, Fop);
        while (feof(Fop) == 0) {
            l = strlen(line);
            if (line[0] == '#' || l < 3) {
                fgets(line, 256, Fop);
                continue;
            }
            if (l > 55) strncpy(key, &line[54], 4);
            else strcpy(key, " PIC");
            key[4] = 0;
            for (i = 0; i < 4; i++) key[i] = toupper(key[i]);
            if (strcmp(key, " PIC") == 0) {
                lp = (NewPick *) realloc(lp, sizeof (NewPick)*(npick + 1));
                memset(&lp[npick], 0, sizeof (NewPick));
                strncpy(str2, line, 18);
                str2[18] = 0;
                n = sscanf(str2, "%s %s %s", lp[npick].stname, lp[npick].phase, lp[npick].quali);

                strncpy(str2, &line[18], 10);
                str2[10] = 0;
                n = sscanf(str2, "%f %d ", &lp[npick].time, &lp[npick].weight);

                strncpy(str2, &line[28], 26);
                str2[26] = 0;
                n = sscanf(str2, "%f %d %*c%s %d", &lp[npick].periode,
                        &lp[npick].amplitude, lp[npick].amptype, &lp[npick].ampuse);
                n = 0;
                if (strcmp(&lp[npick].amptype[1], "WA") == 0) locP->AmpFilt = 'L';
                if (l > 55) { /* amplitude channel and time: col 60-70 */
                    memset(str2, 0, 256);
                    strncpy(str2, &line[60], 11);
                    if (str2[0] != ' ') {
                        n = sscanf(str2, "%s %f", lp[npick].ampchan, &lp[npick].amptime);
                    } else {
                        n = sscanf(&str2[4], "%f", &lp[npick].amptime);
                    }
                    n = 2;
                }
                if (n < 2) {
                    lp[npick].ampchan[0] = 0;
                    lp[npick].amptime = 0;
                }
                if (l > 71) { /* pick uncertainties: col 72-92 */
                    memset(str2, 0, 256);
                    strcpy(str2, &line[71]);
                    n = sscanf(str2, "%f %f", &lp[npick].tmin, &lp[npick].tmax);
                } else {
                    if (lp[npick].quali[0] == 'I') {
                        lp[npick].tmin = -tmin;
                        lp[npick].tmax = tmin;
                    }
                    if (lp[npick].quali[0] == 'E') {
                        lp[npick].tmin = -tmin * 5;
                        lp[npick].tmax = tmin * 5;
                    }
                    if (lp[npick].quali[0] == 'Q' || lp[npick].quali[0] == ' ') {
                        lp[npick].tmin = -tmin * 10;
                        lp[npick].tmax = tmin * 10;
                    }
                }
                lp[npick].minute = eminute;
                lp[npick].amplitude /= 2;
                for (j = 0; j < nstat; j++) {
                    if (strcmp(station[j].name, lp[npick].stname) == 0) {
                        n = 0;
                        for (i = 0; i < station[j].nsig; i++) {
                            if (station[j].shead[i].s_comp[2] == 'Z') {
                                n = i;
                                if (lp[npick].ampchan[0] == 0) strcpy(lp[npick].ampchan, station[j].shead[i].s_comp);
                            }
                        }
                        lp[npick].dmy = station[j].trace_nr[n] + 1;
                        break;
                    }
                }
                npick++;
            }
            if (strcmp(key, "LOCA") == 0) {
                for (i = 0; i < 20; i++) {
                    if (line[i] == '/' || line[i] == ':') line[i] = ' ';
                }
                n = sscanf(line, "%d %d %d %d %d %*s %d", &yr, &mo, &dy, &hr, &mn, &evnr);
                eminute = juliam4(&yr, &mo, &dy, &hr, &mn);
                strcpy(locP->type, "L");
                tmin = 0.05;
            }
            if (strcmp(key, "REGI") == 0) {
                for (i = 0; i < 20; i++) {
                    if (line[i] == '/' || line[i] == ':') line[i] = ' ';
                }
                n = sscanf(line, "%d %d %d %d %d %*s %d", &yr, &mo, &dy, &hr, &mn, &evnr);
                eminute = juliam4(&yr, &mo, &dy, &hr, &mn);
                strcpy(locP->type, "R");
            }
            if (strcmp(key, "TELE") == 0) {
                for (i = 0; i < 20; i++) {
                    if (line[i] == '/' || line[i] == ':') line[i] = ' ';
                }
                n = sscanf(line, "%d %d %d %d %d %*s %d", &yr, &mo, &dy, &hr, &mn, &evnr);
                eminute = juliam4(&yr, &mo, &dy, &hr, &mn);
                strcpy(locP->type, "T");
            }
            if (strcmp(key, "UNKN") == 0) {
                for (i = 0; i < 20; i++) {
                    if (line[i] == '/' || line[i] == ':') line[i] = ' ';
                }
                n = sscanf(line, "%d %d %d %d %d %*s %d", &yr, &mo, &dy, &hr, &mn, &evnr);
                eminute = juliam4(&yr, &mo, &dy, &hr, &mn);
            }
            if (strcmp(key, "TRIA") == 0) {
                line[54] = 0;
                n = sscanf(line, "%f %f %f %d %d %d %d %d %f %d",
                        &locP->lat, &locP->lon, &locP->depth,
                        &locP->yr, &locP->mo, &locP->dy, &locP->hr, &locP->mn,
                        &locP->sec, &locP->fix);
            }
            if (strcmp(key, "INST") == 0) {
                line[54] = 0;
                n = sscanf(line, "%d %d %f %f %f %d", &locP->st_list, &locP->model,
                        &locP->xnear,
                        &locP->xfar, &locP->pos, &locP->swt);
            }
            fgets(line, 256, Fop);
        }
        fclose(Fop);
        *p = lp;
        *Np = npick;
        NrEvent = evnr;
        return ( 0);
    }
    return ( 1);
}

int ReadFiltersfromFile(file, typ, limit, f, Nf)
char *file, *typ;
int limit;
NewFilter **f;
int *Nf;
{
    FILE *Fop;
    char ftyp[8], zerophase;
    int n, order, nn;
    float low, high;

    NewFilter *lf = NULL;

    nn = 0;
    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    if (Fop != NULL) {
        fgets(line, 256, Fop);
        while (feof(Fop) == 0) {
            if (line[0] != '#') {
                n = sscanf(line, "%s %d %f %f %c", ftyp, &order, &low, &high, &zerophase);
                if (n > 2) {
                    lf = (NewFilter *) realloc(lf, sizeof (NewFilter)*(nn + 1));
                    memset(lf, 0, sizeof (NewFilter));
                    strcpy(lf[nn].typ, ftyp);
                    lf[nn].order = order;
                    lf[nn].low = low;
                    if (n > 3) lf[nn].high = high;
                    else lf[nn].high = 0;
                    if (n > 4) {
                        if (zerophase == 'T' || zerophase == 't') lf[nn].zerophase = True;
                        else lf[nn].zerophase = False;
                    } else lf[nn].zerophase = False;
                    nn++;
                }
            }
            fgets(line, 256, Fop);
        }
        fclose(Fop);
        *Nf = nn;
        *f = lf;
        return ( 0);
    }
    return ( 1);
}

int ReadLocationfromFile(file, typ, limit, l)
char *file, *typ;
int limit;
NewLoc *l;
{
    FILE *Fop;

    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    if (Fop != NULL) {
        fgets(line, 256, Fop);
        while (feof(Fop) == 0) {
            if (line[0] != '#') {
                if (strlen(line) > 2) {
                    strcpy((char *) &(l->locsum), line);
                }
            }
            fgets(line, 256, Fop);
        }
        fclose(Fop);
        return ( 0);
    }
    return ( 1);
}

int ReadCalibrationfromFile(file, typ, limit)
char *file, *typ;
int limit;
{
    FILE *Fop;

    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    if (Fop != NULL) {

        fclose(Fop);
        return ( 0);
    }
    return ( 1);
}

int ReadWIDfromFile(file, typ, limit)
char *file;
char *typ;
int limit;
{
    FILE *Fop;
    int i, j, k, l, ntrace;
    float dtime;

    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    ntrace = 0;
    s_freq_max = 0;
    if (Fop != NULL) {
        fgets(line, 256, Fop);
        while (feof(Fop) == 0) {
            l = strlen(line);
            if (l > 5) {
                if (strncmp(line, "WID2 ", 5) == 0) {
                    l = ReadWIDoneChannel(Fop, line, ntrace);
                    fgets(line, 256, Fop);
                    if (l == 0) ntrace++;
                }
            }
        }
        fclose(Fop);
        /* determine overall start and end time */
        time_maximum = 0;
        sec_maximum = 0;
        time_minimum = 2147483647; /* 2^31 -1 */
        sec_minimum = 0;
        for (i = 0; i < nstat; i++) {
            for (j = 0; j < station[i].nsig; j++) {
                station[i].comp_tab[j] = AddComp2Table(station[i].shead[j].s_comp);
                station[i].available[j] = 1 << station[i].comp_tab[j];
                dtime = (float) (time_minimum - station[i].shead[j].s_start_min) *60
                        + sec_minimum - station[i].shead[j].s_start_sec;
                if (dtime > 0) {
                    time_minimum = station[i].shead[j].s_start_min;
                    sec_minimum = station[i].shead[j].s_start_sec;
                }
                dtime = (float) (time_maximum - station[i].shead[j].s_end_min) *60
                        + sec_maximum - station[i].shead[j].s_end_sec;
                if (dtime < 0) {
                    time_maximum = station[i].shead[j].s_end_min;
                    sec_maximum = station[i].shead[j].s_end_sec;
                }
                if (s_freq_max < station[i].shead[j].s_freq) s_freq_max = station[i].shead[j].s_freq;
                if (station[i].shead[j].s_resol > 0) {
                    l = 1;
                    for (k = 0; k < station[i].shead[j].s_resol - 1; k++) l *= 2;
                    station[i].shead[j].s_saturation = (float) l;
                    if (station[i].shead[j].s_resol < 17) station[i].shead[j].s_saturation *= 0.75;
                }
            }
        }
        return ( 0);
    }
    return ( 1);
}

int ReadWIDoneChannel(fp, in, ntrace)
FILE *fp;
char *in;
int ntrace;
{
    int n, i, l, stNr, coNr, ctNr;
    int yr, mo, dy, hr, mn, sampl, jmin, nbyte;
    float sec, srate, nmvolt, calper, hang, vang, dur;
    char sname[8], scomp[8], saux[8], instr[8], compr[8], netw[20];
    char key[6], *wid;

    /* decode WID2 line */
    for (i = 5; i < 29; i++) {
        if (in[i] == '/' || in[i] == ':') in[i] = ' ';
    }
    n = sscanf(in, "%*s %d %d %d %d %d %f %s %s %s %s %d %f %f %f %s %f %f",
            &yr, &mo, &dy, &hr, &mn, &sec, sname, scomp, saux, compr,
            &sampl, &srate, &nmvolt, &calper, instr, &hang, &vang);
    jmin = juliam4(&yr, &mo, &dy, &hr, &mn);
    stNr = checkStationinList(sname);
    ctNr = checkCompTable(scomp);
    coNr = checkComponent(stNr, scomp);
    memset(&station[stNr].shead[coNr], 0, sizeof (SignalHead));
    station[stNr].comp_tab[coNr] = ctNr;
    station[stNr].trace_nr[coNr] = ntrace;
    strcpy(station[stNr].shead[coNr].s_name, sname);
    strcpy(station[stNr].shead[coNr].s_comp, scomp);
    strcpy(station[stNr].shead[coNr].s_typ, instr);
    strcpy(station[stNr].shead[coNr].s_compress, compr);
    strcpy(station[stNr].shead[coNr].s_recsys, saux);
    if (strcmp(compr, "CM6") == 0) station[stNr].shead[coNr].s_ndiff = 2;
    else station[stNr].shead[coNr].s_ndiff = 0;
    station[stNr].shead[coNr].s_start_min = jmin;
    station[stNr].shead[coNr].s_start_sec = sec;
    dur = (float) sampl / srate;
    sec += dur;
    while (sec > 60) {
        jmin++;
        sec -= 60;
    }
    station[stNr].shead[coNr].s_end_min = jmin;
    station[stNr].shead[coNr].s_end_sec = sec;
    station[stNr].shead[coNr].s_freq = srate;
    station[stNr].shead[coNr].s_nsampl = sampl;
    station[stNr].shead[coNr].s_status = -1;
    if (sampl > MaxSample) MaxSample = sampl;

    /* read WID section */
    wid = (char *) malloc(sizeof (int) * sampl);
    memset(wid, 0, sizeof (int) * sampl);
    nbyte = 0;
    fgets(line, 256, fp);
    strncpy(key, line, 5);
    key[5] = 0;
    if (key[4] == '\n') key[4] = ' ';
    while (1) {
        if (strcmp(key, "DATAS") == 0) {
            l = strlen(line) - 1;
            memcpy(&wid[nbyte], line, l);
            nbyte += l;
            fgets(line, 256, fp);
            if (strncmp(line, "CHK2 ", 5) == 0) strcpy(key, "CHK2 ");
            continue;
        }
        if (strcmp(key, "STA2 ") == 0) {
            n = sscanf(line, "%*s %s %f %f %*s %f %f", netw,
                    &station[stNr].shead[coNr].s_lat, &station[stNr].shead[coNr].s_lon,
                    &station[stNr].shead[coNr].s_elev, &station[stNr].shead[coNr].s_burrial);
            fgets(line, 256, fp);
            strncpy(key, line, 5);
            key[5] = 0;
            if (key[4] == '\n') key[4] = ' ';
            continue;
        }
        if (strcmp(key, "DAT2 ") == 0) {
            fgets(line, 256, fp);
            strcpy(key, "DATAS");
            continue;
        }
        if (strcmp(key, "CHK2 ") == 0) {
            n = sscanf(line, "%*s %ld", &station[stNr].shead[coNr].s_checksm);
            break;
        }
    }
    wid = (char *) realloc(wid, nbyte);
    station[stNr].shead[coNr].s_n_compress = nbyte;
    station[stNr].shead[coNr].s_nblocks = (nbyte + 511) / 512;
    station[stNr].wave[coNr] = wid;
    station[stNr].shead[coNr].s_nbytes = 4;
    if (nbyte > 0) return (0);
    else return (1);
}

int checkCompTable(comp)
char *comp;
{
    int i, nc;
    nc = -1;
    for (i = 0; i < ncomp_table; i++) {
        if (strcmp(comp_table[i], comp) == 0) {
            nc = i;
            break;
        }
    }
    if (nc < 0) {
        strcpy(comp_table[ncomp_table], comp);
        nc = ncomp_table;
        ncomp_table++;
    }
    return ( nc);
}

int checkComponent(nrs, comp)
int nrs;
char *comp;
{
    int i, nn;
    nn = -1;
    for (i = 0; i < station[nrs].nsig; i++) {
        if (strcmp(station[nrs].chanid[i], comp) == 0) {
            nn = i;
            break;
        }
    }
    if (nn < 0) {
        nn = station[nrs].nsig;
        strcpy(station[nrs].chanid[nn], comp);
        station[nrs].nsig++;
    }
    return ( nn);
}

int checkStationinList(sname)
char *sname;
{
    int i, nn;
    nn = -1;
    for (i = 0; i < nstat; i++) {
        if (strcmp(station[i].name, sname) == 0) {
            nn = i;
            break;
        }
    }
    if (nn < 0) {
        nn = nstat;
        strcpy(station[nn].name, sname);
        nstat++;
    }
    return ( nn);
}

int ReadCALfromFile(file, typ, limit)
char *file, *typ;
int limit;
{
    FILE *Fop;
    PAZ_str *lpaz = NULL;
    //FIR_str *fir=NULL;
    Digitizer *digit = NULL;
    char ss[30], tt[30];
    int i, k, n, nc = 0;
    int yr, mo, dy, hr, mn, stage;
    float *poles, *zeros;

    strcpy(str1, file);
    Fop = OpenLimiter(str1, typ, "r", limit);
    if (Fop != NULL) {
        fgets(line, 256, Fop);
        while (feof(Fop) == 0) {
            if (line[0] != '#') {
                if (strlen(line) > 2) {
                    if (strncmp(line, "CAL2 ", 5) == 0) {
                        calib = (CalEntry *) realloc(calib, sizeof (CalEntry)*(ncalib + 1));
                        memset(&calib[ncalib], 0, sizeof (CalEntry));
                        stage = 0;
                        n = sscanf(line, "%*s %s %s", calib[ncalib].name, calib[ncalib].comp);
                        strncpy(ss, &line[15], 4);
                        ss[4] = 0;
                        n = sscanf(ss, "%s", calib[ncalib].auxid);
                        strncpy(ss, &line[28], 22);
                        ss[22] = 0;
                        n = sscanf(ss, "%f %f", &calib[ncalib].nmct, &calib[ncalib].calper);
                        sscanf(&line[53], "%f", &calib[ncalib].srate);
                        strncpy(ss, &line[63], 17);
                        ss[17] = 0;
                        for (i = 0; i < 17; i++) {
                            if (ss[i] == '/' || ss[i] == ':') ss[i] = ' ';
                        }
                        n = sscanf(ss, "%d %d %d %d %d", &yr, &mo, &dy, &hr, &mn);
                        calib[ncalib].time1 = juliam4(&yr, &mo, &dy, &hr, &mn);
                        if (strlen(line) > 82) {
                            strncpy(ss, &line[80], 17);
                            ss[17] = 0;
                            for (i = 0; i < 17; i++) {
                                if (ss[i] == '/' || ss[i] == ':') ss[i] = ' ';
                            }
                            n = sscanf(ss, "%d %d %d %d %d", &yr, &mo, &dy, &hr, &mn);
                        } else {
                            yr += 100;
                        }
                        calib[ncalib].time2 = juliam4(&yr, &mo, &dy, &hr, &mn);
                        nc = ncalib;
                        ncalib++;
                    }
                    if (strncmp(line, "PAZ2 ", 5) == 0) {
                        lpaz = (PAZ_str *) malloc(sizeof (PAZ_str));
                        memset(lpaz, 0, sizeof (PAZ_str));
                        n = sscanf(line, "%*s %d %s %f %d %d",
                                &stage, lpaz->units, &(lpaz->factor), &(lpaz->npoles), &(lpaz->nzeros));
                                // 20100617 AJL  &stage, &(lpaz->units), &(lpaz->factor), &(lpaz->npoles), &(lpaz->nzeros));
                        poles = (float *) malloc(sizeof (float) * lpaz->npoles * 2);
                        memset(poles, 0, sizeof (float) * lpaz->npoles * 2);
                        zeros = (float *) malloc(sizeof (float) * lpaz->nzeros * 2);
                        memset(zeros, 0, sizeof (float) * lpaz->nzeros * 2);
                        for (i = 0, k = 0; i < lpaz->npoles; i++, k += 2) {
                            fgets(line, 256, Fop);
                            sscanf(line, "%f %f", &poles[k], &poles[k + 1]);
                        }
                        for (i = 0, k = 0; i < lpaz->nzeros; i++, k += 2) {
                            fgets(line, 256, Fop);
                            sscanf(line, "%f %f", &zeros[k], &zeros[k + 1]);
                        }
                        lpaz->poles = poles;
                        lpaz->zeros = zeros;
                        calib[nc].stage[stage - 1].paz = lpaz;
                        calib[nc].stage[stage - 1].fir = NULL;
                        calib[nc].stage[stage - 1].digit = NULL;
                        if (stage > calib[nc].Nstage) calib[nc].Nstage = stage;
                    }
                    if (strncmp(line, "FIR2 ", 5) == 0) {
                    }
                    if (strncmp(line, "DIG2 ", 5) == 0) {
                        digit = (Digitizer *) malloc(sizeof (Digitizer));
                        memset(digit, 0, sizeof (Digitizer));
                        n = sscanf(line, "%*s %d %f %f", &stage, &(digit->counts_volt), &(digit->srate));
                        calib[nc].stage[stage - 1].paz = NULL;
                        calib[nc].stage[stage - 1].fir = NULL;
                        calib[nc].stage[stage - 1].digit = digit;
                        if (stage > calib[nc].Nstage) calib[nc].Nstage = stage;
                    }
                    if (strncmp(line, "( Res", 5) == 0) {
                        n = sscanf(line, "%*s %*s %d", &calib[nc].resolution);
                    }
                    if (strncmp(line, "( Vol", 5) == 0) {
                        n = sscanf(line, "%*s %*s %f", &calib[nc].volts);
                    }
                    if (strncmp(line, "( Sco", 5) == 0) {
                        n = sscanf(line, "%*s %*s %f", &calib[nc].sconst);
                    }
                    if (strncmp(line, "( Per", 5) == 0) {
                        n = sscanf(line, "%*s %*s %f", &calib[nc].periode);
                    }
                    if (strncmp(line, "( Dam", 5) == 0) {
                        n = sscanf(line, "%*s %*s %f", &calib[nc].damping);
                    }
                    if (strncmp(line, "( Gai", 5) == 0) {
                        n = sscanf(line, "%*s %*s %f", &calib[nc].gain);
                    }
                    if (strncmp(line, "( Rec", 5) == 0) {
                        n = sscanf(line, "%*s %*s %s", tt);
                        strcpy(calib[nc].recsys, tt);
                    }

                }
            }
            fgets(line, 256, Fop);
        }

        fclose(Fop);
        return ( 0);
    }
    return ( 1);
}

char *Comp6toBuf(in4, nrofint, nrofchar, cbuf)
long *in4;
long nrofint;
long *nrofchar;
char *cbuf;
{
    /* Adjust this statement for your memory model or operating system */

#ifdef UNIX
#define MBS 128*1024*1024
#endif
#ifdef MSDOS
#define MBS 65530
#endif

    static char LookUp[65] ={' ', '+', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
    long n4 = 16;
    long n5 = 32;
    long n5m1 = 31;
    //long n6  = 64;
    long n9 = 512;
    long n10 = 1024;
    long n10m1 = 1023;
    long n14 = 16384;
    long n15 = 32768;
    long n15m1 = 32767;
    long n19 = 524288;
    long n20 = 1048576;
    long n20m1 = 1048575;
    long n24 = 16777216;
    long n25 = 33554432;
    long n25m1 = 33554431;
    //long n28  = 134217728;
    long n28m1 = 134217727;
    long mflag = 32;
    long nflag;
    long a;
    long count = 0;
    long jn;
    long j;
    long nblank;
    int i;

    cbuf = (char *) malloc(1);
    if (cbuf == NULL) {
        i = fprintf(stderr, "Memory allocation error in Comp6toBuf (1) !\n");
        exit(1);
    }
    for (a = 0; a < nrofint; a++) {
        nflag = 1;
        jn = in4[a];
        if (jn < 0) {
            nflag = nflag + n4;
            jn = -jn;
        }
        if (jn >= n4) {
            /* more than 1 byte */
            if (jn >= n9) {
                /* more than 2 bytes */
                if (jn >= n14) {
                    /* more than 3 bytes */
                    if (jn >= n19) {
                        /* more than 4 bytes */
                        if (jn >= n24) {
                            /* more than 5 bytes */
                            if (jn >= n28m1)
                                jn = n28m1;
                            j = jn / n25 + nflag + mflag;
                            if ((count + 1) > MBS) {
                                i = fprintf(stderr, "\nCompression buffer became too large!!\n");
                                exit(1);
                            }
                            if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
                                i = fprintf(stderr, "\nMemory allocation error  in Comp6toBuf (2) !\n");
                                exit(1);
                            }
                            cbuf[count] = LookUp[j];
                            count++;
                            jn = jn & n25m1;
                            nflag = 1;
                        }
                        j = (jn / n20 + nflag + mflag);
                        if ((count + 1) > MBS) {
                            i = fprintf(stderr, "\nCompression buffer became too large!!\n");
                            exit(1);
                        }
                        if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
                            i = fprintf(stderr, "\nMemory allocation error in Comp6toBuf (3) !\n");
                            exit(1);
                        }
                        cbuf[count] = LookUp[j];
                        count++;
                        jn = jn & n20m1;
                        nflag = 1;
                    }
                    j = (jn / n15 + nflag + mflag);
                    if ((count + 1) > MBS) {
                        i = fprintf(stderr, "\nCompression buffer became too large!!\n");
                        exit(1);
                    }
                    if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
                        i = fprintf(stderr, "\nMemory allocation error in Comp6toBuf (4) !\n");
                        exit(1);
                    }
                    cbuf[count] = LookUp[j];
                    count++;
                    jn = jn & n15m1;
                    nflag = 1;
                }
                j = (jn / n10 + nflag + mflag);
                if ((count + 1) > MBS) {
                    i = fprintf(stderr, "\nCompression buffer became too large!!\n");
                    exit(1);
                }
                if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
                    i = fprintf(stderr, "\nMemory allocation error in Comp6toBuf (5) !\n");
                    exit(1);
                }
                cbuf[count] = LookUp[j];
                count++;
                jn = jn & n10m1;
                nflag = 1;
            }
            j = (jn / n5 + nflag + mflag);
            if ((count + 1) > MBS) {
                i = fprintf(stderr, "\nCompression buffer became too large!!\n");
                exit(1);
            }
            if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
                i = fprintf(stderr, "\nMemory allocation error in Comp6toBuf (6) !!\n");
                exit(1);
            }
            cbuf[count] = LookUp[j];
            count++;
            jn = jn & n5m1;
            nflag = 1;
        }
        j = (jn + nflag);
        if ((count + 1) > MBS) {
            i = fprintf(stderr, "\nCompression buffer became too large!!\n");
            exit(1);
        }
        if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
            i = fprintf(stderr, "\nMemory allocation error in Comp6toBuf (7) !\n");
            exit(1);
        }
        cbuf[count] = LookUp[j];
        count++;
    }
    nblank = 80 - (count % 80);
    if (nblank < 2) nblank += 80;
    for (a = 0; a < nblank; a++) {
        if ((count + 1) > MBS) {
            i = fprintf(stderr, "\nCompression buffer became too large!!\n");
            exit(1);
        }
        if ((cbuf = (char *) realloc(cbuf, (size_t) (count + 1))) == NULL) {
            i = fprintf(stderr, "\nMemory allocation error in Comp6toBuf (8) !\n");
            exit(1);
        }
        cbuf[count] = ' ';
        count++;
    }
    (*nrofchar) = count;
    return (cbuf);
}

void Decomp6toArray(lb, lout, iout, cbuf)
long lb;
long *lout;
long *iout;
char *cbuf;
{
    static long ichar[129] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2,
        3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0, 0, 0,
        12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
        28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 0, 0, 0, 0, 0, 0,
        38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
        54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 0, 0, 0, 0, 0, 0};
    static char test[5] = {'1', '2', '3', '4', '\0'};
    char achar[5];
    char aspace = 32;
    char lfeed = 10;
    char cretn = 13;
    long *itest;
    long *inn;
    long imax;
    long isign;
    long jsign;
    long ioflow;
    long joflow;
    long mask1;
    long mask2;
    long ibyte;
    long icount;
    long i;
    long j;
    long k;
    long icheck;
    long itemp;

    imax = (*lout);
    isign = 16;
    ioflow = 32;
    mask1 = 15;
    mask2 = 31;
    ibyte = 3;
    itest = (long *) test;
    inn = (long *) achar;
    /* Work out which way bytes are stored in computer */
    icheck = (long) 256 * (long) 256 * (long) 256 * (long) 52;
    icheck += (long) 256 * (long) 256 * (long) 51;
    icheck += (long) 256 * (long) 50;
    icheck += (long) 49;
    if ((*itest) == icheck) {
        ibyte = 0;
    }
    icount = 0;
    i = 0;
    j = 0;
    achar[4] = '\0';
label1:
    i++;
    achar[(unsigned int) ibyte] = cbuf[(unsigned int) i - 1];
    /* If a carriage return or line feed ignore, if a space then end of data */
    if (achar[(unsigned int) ibyte] == cretn) goto label1;
    if (achar[(unsigned int) ibyte] == lfeed) goto label1;
    if (achar[(unsigned int) ibyte] == aspace) goto label5;
    icount++;
    /* Strip off any higher order byte   */
    k = (*inn) & (long) 127;
    /* Get number representation of input character  */
    (*inn) = ichar[(unsigned int) k];
    /* Get sign bit */
    jsign = (*inn) & isign;
    /* get continuation bit (if any)   */
    joflow = (*inn) & ioflow;
    /* Remove bits we dont want    */
    itemp = (*inn) & mask1;
label2:
    if (joflow == 0) goto label4;
    /* There is another byte in this sample   */
    itemp *= 32;
label3:
    i++;
    achar[(unsigned int) ibyte] = cbuf[(unsigned int) i - 1];
    if (achar[(unsigned int) ibyte] == cretn) goto label3;
    if (achar[(unsigned int) ibyte] == lfeed) goto label3;
    icount++;
    /* Strip off any higher order bits   */
    k = (*inn) & (long) 127;
    (*inn) = ichar[(unsigned int) k];
    /* get continuation bit (if any)   */
    joflow = (*inn) & ioflow;
    k = (*inn) & mask2;
    itemp += k;
    goto label2;
label4:
    if (jsign != 0) itemp *= -1;
    j++;
    if (j > imax) goto label5;
    iout[(unsigned int) j - 1] = itemp;


    if (icount < lb) goto label1;
label5:
    if (j > imax) {
        (*lout) = imax;
    }
}

void MakeDiff(in4, nelements, ndiff)
long *in4;
long nelements;
int ndiff;
{
    int a, b;

    b = 0;
    while (b < ndiff) {
        b++;
        for (a = nelements - 1; a >= 1; a--) {
            in4[a] = in4[a] - in4[a - 1];
        }
    }
}

void RemoveDiff(in4, nelements, ndiff)
long *in4;
long nelements;
int ndiff;
{
    int a, b;

    b = 0;
    while (b < ndiff) {
        b++;
        for (a = 1; a < nelements; a++) {
            in4[a] += in4[a - 1];
        }
    }
}

int WritePicksToFile(file, ext, limit, pick, npick, locp, extp)
char *file;
char *ext;
int limit;
NewPick *pick;
int npick;
NewLocPar locp;
int extp; /* =0: short pick line; =1 long pick line */
{
    FILE *fp = NULL;
    char type[10], *lfix_mag;
    int yr, mo, dy, hr, mn, jmin;
    int i, k, m, l;
    float sec, tmin = 0.0;

    if (strlen(file) == 0) {
        fp = stdout;
    } else {
        if ((fp = OpenLimiter(file, ext, "w", limit)) == NULL) {
            l = fprintf(stderr, "unable to open file: %s.%s\n", file, ext);
            return (1);
        }
    }

    k = sortPicks(pick, npick);

    /* INST */
    l = sprintf(line, "%2d%2d%8.2f%8.2f%8.3f%2d%28s",
            locp.model, locp.st_list, locp.xnear,
            locp.xfar, locp.pos, locp.swt, "INST");
    l = fprintf(fp, "%s\n", line);
    /* TRIAL */
    l = sprintf(line, "%8.3f%10.3f%7.3f %04d %02d %02d %02d %02d %5.2f 0%9s",
            locp.lat, locp.lon, locp.depth, locp.yr, locp.mo, locp.dy, locp.hr, locp.mn, locp.sec, "TRIAL");
    l = fprintf(fp, "%s\n", line);
    if (strlen(author) > 0) {
        fprintf(fp, "%-54sAUTOR\n", author);
    }
    if (fabs(locp.dray) > 0.00001 || fabs(locp.daz) > 0.00001) {
        sprintf(line, "%7.4f %7.3f%39cCORR\n", locp.dray, locp.daz, ' ');
        fputs(line, fp);
    }
    lfix_mag = FixMagnitude(line);
    if (strlen(lfix_mag) > 0) i = fprintf(fp, "%s\n", lfix_mag);
    /* Tele Local Regional*/
    if (locp.type[0] == 'T') {
        strcpy(type, "Tele");
        tmin = 0.2;
    } else if (locp.type[0] == 'R') {
        strcpy(type, "Regi");
        tmin = 0.2;
    } else if (locp.type[0] == 'L') {
        strcpy(type, "Loca");
        tmin = 0.05;
    } else {
        strcpy(type, "Unkn");
    }
    jmin = pick[0].minute;
    datum4(&jmin, &yr, &mo, &dy, &hr, &mn);
    l = sprintf(line, "%04d/%02d/%02d %02d:%02d%42s", yr, mo, dy, hr, mn, type);
    if (NrEvent > 0) {
        for (i = strlen(line); i < 78; i++) line[i] = ' ';
        sprintf(&line[69], "%9d", NrEvent);
    }
    l = fprintf(fp, "%s\n", line);
    /* triggers/picks */
    for (i = 0; i < npick; i++) {
        m = sort_idx[i];
        sec = (pick[m].minute - jmin)*60.0 + pick[m].time;
        if (pick[m].phase[0] == 0) strcpy(pick[m].phase, " ");
        if (pick[m].quali[0] == 0) strcpy(pick[m].quali, " ");
        if (pick[m].ampchan[0] == 0) strcpy(pick[m].ampchan, " ");
        if (pick[m].recsys[0] == 0) strcpy(pick[m].recsys, " ");
        if (pick[m].amptype[0] == 0) strcpy(pick[m].amptype, &pick[m].recsys[1]);
        if (fabs(pick[m].tmin) < 0.00001 && fabs(pick[m].tmax) < 0.00001) {
            if (pick[m].quali[0] == 'I') {
                pick[m].tmin = -tmin;
                pick[m].tmax = tmin;
            }
            if (pick[m].quali[0] == 'E') {
                pick[m].tmin = -tmin * 5;
                pick[m].tmax = tmin * 5;
            }
            if (pick[m].quali[0] == 'Q' || pick[m].quali[0] == ' ') {
                pick[m].tmin = -tmin * 10;
                pick[m].tmax = tmin * 10;
            }
        }
        l = sprintf(line, "%-8s%-8s%-2s %7.3f%2d%5.2f%9d %c%-8s%2d Pick %-3s %7.3f %7.3f %7.3f",
                pick[m].stname,
                pick[m].phase,
                pick[m].quali,
                sec,
                pick[m].weight,
                pick[m].periode,
                pick[m].amplitude * 2,
                pick[m].recsys[0],
                pick[m].amptype,
                pick[m].ampuse,
                pick[m].ampchan,
                pick[m].amptime,
                pick[m].tmin,
                pick[m].tmax);
        if (extp == 0) {
            line[54] = 0;
        }

        l = fprintf(fp, "%s\n", line);
    }
    l = fprintf(fp, "%58s\n", "SKIP");
    fclose(fp);
    return (0);
}

char *FixMagnitude(line)
char *line;
{
    char mag_string[30];
    //char str_fixmag[80];
    int l;

    line[0] = 0;
    mag_string[0] = 0;
    if (mag_fixed.fix_ml > 0.0) {
        l = sprintf(mag_string, "Ml=%.2f ", mag_fixed.fix_ml);
        strcat(line, mag_string);
    }
    mag_string[0] = 0;
    if (mag_fixed.fix_mb > 0.0) {
        l = sprintf(mag_string, "mb=%.2f ", mag_fixed.fix_mb);
        strcat(line, mag_string);
    }
    mag_string[0] = 0;
    if (mag_fixed.fix_ms > 0.0) {
        if (strlen(mag_fixed.fix_ms_str) > 0) {
            trim(mag_fixed.fix_ms_str, mag_string);
        }
        strcat(line, mag_string);
        l = sprintf(mag_string, "=%.2f", mag_fixed.fix_ms);
        strcat(line, mag_string);
    }
    if (strlen(line) > 0) {
        expandString(line, 54, ' ');
        strcat(line, "MAGNITUDE");
    }
    return (line);
}

int MagnitudeToFile(file, ext, limit)
char *file;
char *ext;
int limit;
{
    FILE *fp = NULL;
    //char type[10];
    int l, i;

    line[0] = 0;

    if (mag_fixed.fix_ml != 0) {
        l = strlen(line);
        i = sprintf(&line[l], "Ml=%.2f ", mag_fixed.fix_ml);
    }
    if (mag_fixed.fix_mb != 0) {
        l = strlen(line);
        i = sprintf(&line[l], "mb=%.2f ", mag_fixed.fix_mb);
    }
    if (mag_fixed.fix_ms != 0) {
        l = strlen(line);
        i = sprintf(&line[l], "%s=%.2f ", mag_fixed.fix_ms_str, mag_fixed.fix_ms);
    }
    if (strlen(line) > 0) {
        expandString(line, 54, ' ');
        strcat(line, "MAGNITUDE");
        if (strlen(file) == 0) {
            fp = stdout;
        } else {
            if ((fp = OpenLimiter(file, ext, "w", limit)) == NULL) {
                i = fprintf(stderr, "unable to open file: %s.%s \n", file, ext);
                return (1);
            }
        }
        i = fprintf(fp, "%s\n", line);
        fclose(fp);
    }
    return (0);
}

int UpdateLocationOutput(file, list, type)
char *file;
char *list;
char *type;
{
    //char what[10], extens[40], nfile[256];
    char what[10], nfile[256];
    int i, loctype;
    FILE *fpin, *fpout;

    strcpy(what, type);
    for (i = 0; i < strlen(what); i++) what[i] = toupper(what[i]);
    if (strlen(file) <= 0) strcpy(nfile, cleanName);
    else strcpy(nfile, file);
    i = -1;
    if (IsGSEFormat == True) {
        if (what[0] == 'M') loctype = MPDE;
        else if (what[0] == 'A') loctype = APDE;
        else loctype = NOEXT;
        if ((fpout = OpenLimiter(nfile, fext[loctype].extension, "w", fext[loctype].limit)) == NULL) {
            i = fprintf(stderr, "UpdateLocationOutput: unable to open file: %s.%s \n", nfile, fext[loctype].extension);
            return (1);
        }
        if ((fpin = fopen(list, "r")) == NULL) {
            i = fprintf(stderr, "UpdateLocationOutput: unable to open file: %s\n", list);
            return (2);
        }
        fgets(line, 256, fpin);
        while (feof(fpin) == 0) {
            fprintf(fpout, "%s", line);
            fgets(line, 256, fpin);
        }
        fclose(fpin);
        fclose(fpout);
    }
    return (i);
}

int UpdateLocation(file, loc, type)
char *file;
NewLoc *loc;
char *type;
{
    char what[10];
    int i, loctype;

    strcpy(what, type);
    for (i = 0; i < strlen(what); i++) what[i] = toupper(what[i]);
    if (strlen(file) <= 0) strcpy(str1, cleanName);
    else strcpy(str1, file);
    i = -1;
    if (IsGSEFormat == True) {
        if (what[0] == 'M') loctype = MLOC;
        else if (what[0] == 'A') loctype = ALOC;
        else loctype = NOEXT;
        i = WriteLocToFile(str1, fext[loctype].extension, fext[loctype].limit, loc->locsum, "w");
        i = MagnitudeToFile(str1, fext[MAGNITUDE].extension, fext[MAGNITUDE].limit);
    } else {
        i = SaveLocKPFile(str1, what, loc);
    }
    return (i);
}

int SaveLocKPFile(file, what, loc)
char *file;
char *what;
NewLoc *loc;
{
    FILE *fp_data;
    int i, k, nbytes;
    Master m_head;

    if (strlen(what) == 0) return (1);
    if ((fp_data = fopen(file, "rb+")) == NULL) {
        i = fprintf(stderr, "unable to open file: %s\n", file);
        return (1);
    }
    /****************************************
     * read master record:
     ****************************************/
    nbytes = fread(&m_head, sizeof (Master), 1, fp_data);
    if (m_head.loc_rec != 0) {
        i = fseek(fp_data, (long) (m_head.loc_rec - 1)*512, SEEK_SET);
        nbytes = fread(&location, sizeof (kpLocRec), 1, fp_data);
    } else {
        m_head.loc_rec = m_head.lastrec + 1;
        m_head.lastrec++;
    }
    if (what[0] == 'A') {
        memset(location.auto_loc, 32, 80);
        strncpy(location.auto_loc, loc->locsum, strlen(loc->locsum));
        memset(location.auto_reg, 32, 80);
        strncpy(location.auto_reg, loc->region, strlen(loc->region));
    } else if (what[0] == 'M') {
        memset(location.manu_loc, 32, 80);
        strncpy(location.manu_loc, loc->locsum, strlen(loc->locsum));
        memset(location.manu_reg, 32, 80);
        strncpy(location.manu_reg, loc->region, strlen(loc->region));
    } else return (2);

    location.fix_ml = mag_fixed.fix_ml;
    location.fix_mb = mag_fixed.fix_mb;
    location.fix_ms = mag_fixed.fix_ms;
    strcpy(location.fix_ms_str, mag_fixed.fix_ms_str);
    i = fseek(fp_data, (long) (m_head.loc_rec - 1)*512, SEEK_SET);
    i = fwrite(&location, sizeof (kpLocRec), 1, fp_data);
    i = fseek(fp_data, 0L, SEEK_SET);
    if (NrEvent > 10) {
        i = NrEvent / 10000;
        k = NrEvent - i * 10000;
        m_head.ev_nr = (unsigned short) k;
    }
    i = fwrite(&m_head, sizeof (Master), 1, fp_data);
    fclose(fp_data);
    return (0);
}

int WriteLocToFile(file, ext, limit, loc, mode)
char *file;
char *ext;
int limit;
char *loc;
char *mode;
{
    FILE *fp = NULL;
    char type[10];
    int i;

    if (loc[0] != 0) {
        if (strlen(file) == 0) {
            fp = stdout;
        } else {
            strcpy(type, mode);
            if (strlen(type) <= 0) strcpy(type, "a");
            if ((fp = OpenLimiter(file, ext, type, limit)) == NULL) {
                i = fprintf(stderr, "unable to open file: %s.%s with mode <%s>\n", file, ext, type);
                return (1);
            }
        }
        i = fprintf(fp, "%s\n", loc);
        fclose(fp);
    }
    return (0);
}

int WriteFilterToFile(file, ext, limit, filter, nf)
char *file;
char *ext;
int limit;
NewFilter *filter;
int nf;
{
    int i, l;
    char zero;
    FILE *fp;

    if ((fp = OpenLimiter(file, ext, "w", limit)) == NULL) {
        l = fprintf(stderr, "WriteFilterToFile fopen on %s.%s failed\n", file, ext);
        return (1);
    }
    for (i = 0; i < nf; i++) {
        if (filter[i].zerophase) zero = 'Z';
        else zero = ' ';
        l = fprintf(fp, "%s %d %.3f %.3f %c\n", filter[i].typ, filter[i].order,
                filter[i].low, filter[i].high, zero);
    }
    fclose(fp);
    return (0);
}

int WriteSignalStatusToFile(file, ext, limit)
char *file;
char *ext;
int limit;
{
    int i, j, k, l;
    FILE *fp = NULL;
    float *sz = NULL;
    int Nsz = 0, nsmp;

    for (i = 0; i < nstat; i++) {
        for (j = 0; j < station[i].nsig; j++) {
            if (Nsz < station[i].shead[j].s_nsampl) {
                Nsz = station[i].shead[j].s_nsampl + 100;
                sz = (float *) realloc(sz, sizeof (float) * Nsz);
            }
            l = ConvertToFloat(&station[i].shead[j], station[i].wave[j], sz, &nsmp);
            if (station[i].shead[j].status_out != NULL) {
                for (k = 0; k < station[i].shead[j].Nstatus_out; k++) {
                    if (fp == NULL) {
                        if ((fp = OpenLimiter(file, ext, "w", limit)) == NULL) {
                            l = fprintf(stderr, "WriteSignalStatusToFile fopen on %s.%s failed\n", file, ext);
                            return (1);
                        }
                    }
                    fprintf(fp, "%s %s %d %d\n", station[i].shead[j].s_name,
                            station[i].shead[j].s_comp,
                            station[i].shead[j].status_out[k],
                            station[i].shead[j].status_back[k]);
                }
            }
        }
    }
    if (sz != NULL) free(sz);
    if (fp != NULL) fclose(fp);
    return (0);
}

int ConvertToFloat(shead, wave, si, nsamp)
SignalHead *shead;
char *wave;
float *si;
int *nsamp;
{

    // AJL 20040518
    long *buffer;
    long ns, nn;
    //  int *buffer;
    //   int i, j, k, l, ist, ns, nn;
    int i, j, k, l, ist;
    // AJL
    int signum, kshft, io1;
    char *d_char;
    short *d_short;
    int *d_int;
    double two;
    //double two, exp, resl;

    if (shead->s_volts >= 0.0) signum = 1;
    else signum = -1;
    ns = shead->s_n_compress;
    nn = shead->s_nsampl;
    buffer = NULL;
    // AJL 20040518
    if (ns < nn) {
        buffer = (long *) malloc(sizeof (long) * nn);
    } else {
        buffer = (long *) malloc(sizeof (long) * ns);
    }
    /*
    if ( ns < nn ) {
      buffer= (int *) malloc(sizeof(int) * nn);
    }
    else {
      buffer= (int *) malloc(sizeof(int) * ns);
    }*/
    // AJL
    if (buffer == NULL) return (3);
    if (strncmp(shead->s_compress, "BIN", 3) == 0) { /* no compression */
        switch (shead->s_nbytes) {
            case 1:
            {
                d_char = wave;
                for (k = 0; k < shead->s_nsampl; k++) buffer[k] = d_char[k];
                break;
            }
            case 2:
            {
                d_short = (short *) wave;
                for (k = 0; k < shead->s_nsampl; k++) buffer[k] = d_short[k];
                break;
            }
            case 4:
            {
                d_int = (int *) wave;
                for (k = 0; k < shead->s_nsampl; k++) buffer[k] = d_int[k];
                break;
            }
            default:
            {
                l = fprintf(stderr, "can not input %d  in BIN format\n", shead->s_nbytes);
                *nsamp = 0;
                free(buffer);
                return (2);
            }
        }
    } else {
        if (strncmp(shead->s_compress, "CM6", 4) == 0) { /* 6-bit compression */
#ifdef Matlab
            mexPrintf("%d \n", nn);
            mexPrintf("%10s\n", sb);
            Decomp6toArray(ns, &nn, buffer, wave);
            mexPrintf("%d \n", nn);
            mexPrintf("%d %d %d %d\n", buffer[0], buffer[1], buffer[2], buffer[3]);
#else
            Decomp6toArray(ns, &nn, buffer, wave);
#endif	   
        } else {
            /* CMP7 and CMP8 not supported */
            l = fprintf(stderr, " data storage type %4s not supported!\n", shead->s_compress);
            shead->s_nsampl = 0;
            *nsamp = 0;
            free(buffer);
            return (3);
        }
    }
    /* remove differences */
    if (shead->s_ndiff) {
        RemoveDiff(buffer, nn, shead->s_ndiff);
    }
    ist = 1 << shead->s_status;
    /**
    kshft= 2**(shead->s_position-shead->s_resol + 1)
     **/
    two = 2;
    k = shead->s_position - shead->s_resol + 1;
    if (shead->s_position == 0) {
        kshft = 0;
    } else {
        kshft = 1;
        for (i = 0; i < k; i++) kshft *= 2;
    }
    if (shead->s_status >= 0) {
        if (shead->status_out != NULL) {
            free(shead->status_out);
            shead->status_out = NULL;
            shead->status_back = NULL;
        }
        io1 = 0;
        shead->Nstatus_out = 0;
        for (i = 0; i < nn; i++) {
            k = buffer[i];
            if (k & ist) {
                if (io1 == 0) {
                    shead->status_out = (int *) realloc(shead->status_out, sizeof (int) *(shead->Nstatus_out + 1));
                    shead->status_back = (int *) realloc(shead->status_back, sizeof (int) *(shead->Nstatus_out + 1));
                    shead->status_out[shead->Nstatus_out] = i;
                    io1 = 1;
                }
            } else {
                if (io1 > 0) {
                    shead->status_back[shead->Nstatus_out] = i;
                    io1 = 0;
                    shead->Nstatus_out++;
                }
            }
            buffer[i] = signum * (k / kshft);
        }
        if (io1 > 0) {
            shead->status_back[shead->Nstatus_out] = nn - 1;
            shead->Nstatus_out++;
        }
    } else if (kshft > 1) {
        for (i = 0; i < nn; i++) buffer[i] = signum * (buffer[i] / kshft);
    }
    *nsamp = nn;
    if (DoStatusCheck) {
        for (j = 0; j < shead->Nstatus_out; j++) {
            for (i = shead->status_out[j]; i < shead->status_back[j]; i++) buffer[i] = 0;
        }
    }
    for (i = 0; i < nn; i++) si[i] = (float) buffer[i];
#ifdef Matlab
    mexPrintf("%d %d %d %d\n", buffer[0], buffer[1], buffer[2], buffer[3]);
    mexPrintf("%f %f %f %f\n", (double) si[0], (double) si[1], (double) si[2], (double) si[3]);
#endif
    free(buffer);
    shead->s_nsampl = nn;
    return (0);
}

int ConvertToInt(shead, wave, si, nsamp)
SignalHead *shead;
char *wave;
int *si;
int *nsamp;
{
    // AJL 20040518
    long *buffer;
    long ns, nn;
    //  int *buffer;
    //   int i, j, k, l, ist, ns, nn;
    int i, j, k, l, ist;
    // AJL
    int signum, kshft, io1;
    char *d_char;
    short *d_short;
    int *d_int;
    double two;
    //double two, exp, resl;

    if (shead->s_volts >= 0.0) signum = 1;
    else signum = -1;
    ns = shead->s_n_compress;
    nn = shead->s_nsampl;
    buffer = NULL;
    // AJL 20040518
    if (ns < nn) {
        buffer = (long *) malloc(sizeof (long) * nn);
    } else {
        buffer = (long *) malloc(sizeof (long) * ns);
    }
    /*
    if ( ns < nn ) {
      buffer= (int *) malloc(sizeof(int) * nn);
    }
    else {
      buffer= (int *) malloc(sizeof(int) * ns);
    }*/
    // AJL
    if (buffer == NULL) return (3);
    if (strncmp(shead->s_compress, "BIN", 3) == 0) { /* no compression */
        switch (shead->s_nbytes) {
            case 1:
            {
                d_char = wave;
                for (k = 0; k < shead->s_nsampl; k++) buffer[k] = d_char[k];
                break;
            }
            case 2:
            {
                d_short = (short *) wave;
                for (k = 0; k < shead->s_nsampl; k++) buffer[k] = d_short[k];
                break;
            }
            case 4:
            {
                d_int = (int *) wave;
                for (k = 0; k < shead->s_nsampl; k++) buffer[k] = d_int[k];
                break;
            }
            default:
            {
                l = fprintf(stderr, "can not input %d  in BIN format\n", shead->s_nbytes);
                *nsamp = 0;
                free(buffer);
                return (2);
            }
        }
    } else {
        if (strncmp(shead->s_compress, "CM6", 4) == 0) { /* 6-bit compression */
            Decomp6toArray(ns, &nn, buffer, wave);
        } else {
            /* CMP7 and CMP8 not supported */
            l = fprintf(stderr, " data storage type %4s not supported!\n", shead->s_compress);
            shead->s_nsampl = 0;
            *nsamp = 0;
            free(buffer);
            return (3);
        }
    }
    /* remove differences */
    if (shead->s_ndiff) {
        RemoveDiff(buffer, nn, shead->s_ndiff);
    }
    ist = 1 << shead->s_status;
    /**
    kshft= 2**(shead->s_position-shead->s_resol + 1)
     **/
    two = 2;
    k = shead->s_position - shead->s_resol + 1;
    if (shead->s_position == 0) {
        kshft = 0;
    } else {
        kshft = 1;
        for (i = 0; i < k; i++) kshft *= 2;
    }
    if (shead->s_status >= 0) {
        if (shead->status_out != NULL) {
            free(shead->status_out);
            shead->status_out = NULL;
            shead->status_back = NULL;
        }
        io1 = 0;
        shead->Nstatus_out = 0;
        for (i = 0; i < nn; i++) {
            k = buffer[i];
            if (k & ist) {
                if (io1 == 0) {
                    shead->status_out = (int *) realloc(shead->status_out, sizeof (int) *(shead->Nstatus_out + 1));
                    shead->status_back = (int *) realloc(shead->status_back, sizeof (int) *(shead->Nstatus_out + 1));
                    shead->status_out[shead->Nstatus_out] = i;
                    io1 = 1;
                }
            } else {
                if (io1 > 0) {
                    shead->status_back[shead->Nstatus_out] = i;
                    io1 = 0;
                    shead->Nstatus_out++;
                }
            }
            buffer[i] = signum * (k / kshft);
        }
        if (io1 > 0) {
            shead->status_back[shead->Nstatus_out] = nn - 1;
            shead->Nstatus_out++;
        }
    } else if (kshft > 1) {
        for (i = 0; i < nn; i++) buffer[i] = signum * (buffer[i] / kshft);
    }
    *nsamp = nn;
    if (DoStatusCheck) {
        for (j = 0; j < shead->Nstatus_out; j++) {
            for (i = shead->status_out[j]; i < shead->status_back[j]; i++) buffer[i] = 0;
        }
    }
    memcpy(si, buffer, nn * sizeof (int));
    free(buffer);
    shead->s_nsampl = nn;
    return (0);
}

float AsignalC(signal, n)
float *signal;
int n;
{
    float jm, js;
    int i;

    jm = 0.0;
    for (i = 0; i < n; i++) {
        js = fabs(signal[i]);
        if (js > jm) jm = js;
    }
    if (jm < 1.0) jm = 1.0;
    return (jm);
}

int sortPicks(pick, npick)
NewPick *pick;
int npick;
{
    int i, j, minmin, *lsort_idx, m;

    lsort_idx = NULL;
    if (npick <= 0) {
        fprintf(stderr, "sortPicks with npick= %d\n", npick);
        return (1);
    }
    sort_idx = (int *) realloc(sort_idx, npick * sizeof (int));
    lsort_idx = (int *) realloc(lsort_idx, npick * sizeof (int));
    sort_time = (float *) realloc(sort_time, npick * sizeof (float));

    minmin = 2147483647;
    for (i = 0; i < npick; i++) {
        if (minmin > pick[i].minute) minmin = pick[i].minute;
    }
    for (i = 0; i < npick; i++) {
        sort_time[i] = (pick[i].minute - minmin)*60 + pick[i].time;
        lsort_idx[i] = i;
    }
    bubble_sort(sort_time, lsort_idx, npick);
    m = 0;
    for (i = 0; i < npick; i++) {
        if (lsort_idx[i] < 0) continue;
        sort_idx[m] = lsort_idx[i];
        m++;
        for (j = i + 1; j < npick; j++) {
            if (lsort_idx[j] < 0) continue;
            if (strcmp(pick[lsort_idx[i]].stname, pick[lsort_idx[j]].stname) == 0) {
                sort_idx[m] = lsort_idx[j];
                m++;
                lsort_idx[j] = -1;
            }
        }
    }
    free(lsort_idx);
    return (0);
}

void bubble_sort(sort, index, nsort)
float sort[];
int index[];
int nsort;
{
    int i, j, itemp;
    float temp;

    nsort--;
    for (i = 0; i < nsort; i++) {
        for (j = nsort; i < j; j--) {
            if (sort[j - 1] > sort[j]) {
                temp = sort[j - 1];
                sort[j - 1] = sort[j];
                sort[j] = temp;
                itemp = index[j - 1];
                index[j - 1] = index[j];
                index[j] = itemp;
            }
        }
    }
}

void bubble_sortI(sort, index, nsort)
int sort[];
int index[];
int nsort;
{
    int i, j, itemp;
    int temp;

    nsort--;
    for (i = 0; i < nsort; i++) {
        for (j = nsort; i < j; j--) {
            if (sort[j - 1] > sort[j]) {
                temp = sort[j - 1];
                sort[j - 1] = sort[j];
                sort[j] = temp;
                itemp = index[j - 1];
                index[j - 1] = index[j];
                index[j] = itemp;
            }
        }
    }
}

/************************************
Decode KP Files
 *************************************/

int DecodeKPfile(file_data, file_length)
unsigned char *file_data;
int file_length;
{
    Master *mmrec;
    int yr, mo, day, hr, mn;

    Nmaster = 0;
    mrec = (R256 *) realloc(mrec, sizeof (R256) * (Nmaster + 1));
    if (mrec == NULL) return (1);
    memcpy(mrec[Nmaster], file_data, 512);
    Nmaster++;
    while (mrec[Nmaster - 1][255] != 0) {
        mrec = (R256 *) realloc(mrec, sizeof (R256) * (Nmaster + 1));
        if (mrec == NULL) return (1);
        memcpy(mrec[Nmaster], &file_data[mrec[Nmaster][255]*512 - 512], 512);
        Nmaster++;
    }
    mmrec = (Master *) mrec[0];
    datum4(&mmrec->starttime, &yr, &mo, &day, &hr, &mn);
    /*   Initialize(); */
    NrEvent = mmrec->ev_nr;
    yr *= 10000;
    NrEvent += yr;
    DecodeKPsignalHead(mmrec, mrec, Nmaster, file_data, file_length);
    DecodeKPTrigger(mmrec, file_data, file_length);
    DecodeKPaPicks(mmrec, file_data, file_length);
    DecodeKPmPicks(mmrec, file_data, file_length);
    DecodeKPamLoc(mmrec, file_data, file_length);
    DecodeKPFilter(mmrec, file_data, file_length);
    return ( 0);
}

void DecodeKPTrigger(m_head, file_data, file_length)
Master *m_head;
unsigned char *file_data;
int file_length;
{
    int i, j;
    if (m_head->nmx_tr_rec > 0) {
        NnmxTrigger = ReadNewPicks(m_head->nmx_tr_rec, file_data, file_length, 3, &nmxTrigger, &locNmxT);
        for (i = 0; i < NnmxTrigger; i++) {
            for (j = 0; j < nstat; j++) {
                if (strcmp(station[j].name, nmxTrigger[i].stname) == 0) {
                    nmxTrigger[i].dmy = station[j].trace_nr[2] + 1;
                    break;
                }
            }
        }
    } else if (m_head->trigger_rec > 0) {
    }
}

void old2new_picks(picks, npicks)
NewPick *picks;
int npicks;
{
    int n;

    for (n = 0; n < npicks; n++) {
        switch (picks[n].quali[0]) {
            case 'I':
            {
                switch (picks[n].weight) {
                    case 4:
                    {
                        picks[n].weight = 0;
                    }
                        break;
                    default:
                    {
                        picks[n].weight = 1;
                    }
                        break;
                }
            }
                break;
            default:
            {
                switch (picks[n].weight) {
                    case 0:
                    case 1:
                    {
                        picks[n].quali[0] = 'I';
                        picks[n].weight = 1;
                    }
                        break;
                    case 2:
                    case 3:
                    {
                        picks[n].quali[0] = 'E';
                        picks[n].weight = 1;
                    }
                        break;
                    default:
                    {
                        picks[n].quali[0] = 'Q';
                        picks[n].weight = 0;
                    }
                }
            }
        }
        if (picks[n].phase[0] == 'U') picks[n].weight = 0;
    }
}

int ReadOldPicks(nrec, file_data, file_length, pick_type, npicksp, lpar)
int nrec;
unsigned char *file_data;
int file_length;
unsigned short pick_type;
NewPick **npicksp;
NewLocPar *lpar;
{

    typedef struct {
        short p_chan;
        char p_typ[4];
        char p_dir[2];
        short p_wt;
        short p_arr;
        short p_noise_amp;
        short p_noise_scale;
        short p_noise_time;
        short p_phase_amp;
        short p_phase_scale;
        short p_phase_time;
        short p_max_amp;
        short p_max_scale;
        short p_max_time;
    } one_pick;

    typedef struct {
        short p_picks;
        short p_tupls;
        long p_process;
        short p_d1[11];
        one_pick st_pick[16];
        unsigned short p_next;
    } pick_rec;

    pick_rec *prec;
    NewPick *picks = NULL;
    SignalHead *shead;
    int baddr, js, ks, ktupl, npicks;
    int j, jstat, jst, nsmp;
    float parr;

    npicks = 0;
    ks = 0;
    while (nrec > 0) {
        baddr = (nrec - 1) * 512;
        if (baddr > file_length) {
            break;
        }
        prec = (pick_rec *) & file_data[baddr];
        js = prec->p_picks;
        ktupl = 256 / prec->p_tupls - 1;
        for (j = 0; j < ktupl; j++) {
            if (ks < js) {
                picks = (NewPick *) realloc(picks, sizeof (NewPick)*(npicks + 1));
                memset(&picks[npicks], 0, sizeof (NewPick));
                jstat = prec->st_pick[j].p_chan - 1;
                jst = channel_to_station[jstat];
                nsmp = prec->st_pick[j].p_arr - 1;
                shead = &station[jst].shead[channel_to_compnr[jstat]];

                parr = shead->s_start_sec;
                parr = parr + (float) nsmp / shead->s_freq;

                strcpy(picks[npicks].stname, shead->s_name);
                picks[npicks].time = parr;
                picks[npicks].minute = shead->s_start_min;
                picks[npicks].weight = prec->st_pick[j].p_wt;
                strncpy(picks[npicks].quali, prec->st_pick[j].p_dir, 2);
                strncpy(picks[npicks].phase, prec->st_pick[j].p_typ, 4);
                trimlencc(picks[npicks].phase, 8, 0);
                strcpy(picks[npicks].ampchan, shead->s_comp);
                if (picks[npicks].dmy <= 0) {
                    picks[npicks].dmy = jstat + 1;
                }
                picks[npicks].amptime = parr + (float) (prec->st_pick[j].p_phase_time) / shead->s_freq;
                picks[npicks].amplitude = prec->st_pick[j].p_phase_amp;
                picks[npicks].periode = 0.5;
                npicks++;
                ks++;
            }
        }
        nrec = prec->p_next;
    }
    if ((pick_type & (short) 2) == 0) {
        old2new_picks(picks, npicks);
    }
    DefaultLocPar(lpar);
    lpar->n_picks = npicks;
    *npicksp = picks;
    return (npicks);
}

int ReadNewPicks(nrec, file_data, file_length, pick_type, npicksp, lpar)
int nrec;
unsigned char *file_data;
int file_length;
unsigned short pick_type;
NewPick **npicksp;
NewLocPar *lpar;
{
    int baddr, first, ks, npicks;
    int js = 0;
    int j, jstat;
    float tmin;

    typedef struct {
        short p_chan;
        char p_typ[4];
        char p_dir[2];
        short p_wt;
        short p_arr;
        short p_noise_amp;
        short p_noise_scale;
        short p_noise_time;
        short p_phase_amp;
        short p_phase_scale;
        short p_phase_time;
        short p_max_amp;
        short p_max_scale;
        short p_max_time;
    } one_pick;

    typedef struct {
        short p_picks;
        short p_tupls;
        long p_process;
        short p_d1[11];
        one_pick st_pick[16];
        unsigned short p_next;
    } pick_rec;

    mpick_rec *mprec = NULL;
    //pick_rec   *prec=NULL;
    NewPick *picks = NULL;

    tmin = 0.2;
    first = -1;
    ks = 0;
    npicks = 0;
    while (nrec != 0) {
        baddr = (nrec - 1) * 512;
        mprec = (mpick_rec *) & file_data[baddr];
        if (first < 0) {
            lpar->n_picks = mprec->p[0].par.n_picks;
            lpar->st_list = mprec->p[0].par.st_list - 48;
            lpar->model = mprec->p[0].par.model - 48;
            lpar->lat = mprec->p[0].par.lat;
            lpar->lon = mprec->p[0].par.lon;
            lpar->depth = mprec->p[0].par.depth;
            lpar->xnear = mprec->p[0].par.xnear;
            lpar->xfar = mprec->p[0].par.xfar;
            lpar->pos = mprec->p[0].par.pos;
            lpar->dray = mprec->p[0].par.dray;
            lpar->daz = mprec->p[0].par.daz;
            lpar->fix = mprec->p[0].par.fix;
            lpar->swt = mprec->p[0].par.swt - 48;
            lpar->AmpFilt = mprec->p[0].par.AmpFilt;
            lpar->NoMag = mprec->p[0].par.NoMag;
            lpar->dumy2 = mprec->p[0].par.dumy2;
            strcpy(lpar->type, mprec->p[0].par.type);
            if (lpar->type[0] == 'L' || lpar->type[0] == 'l') tmin = 0.05;
            js = lpar->n_picks;
            first = 1;
        } else {
            first = 0;
        }
        for (j = first; j < 10; j++) {
            if (ks < js) {
                jstat = mprec->p[j].pick.chan - 1;
                if (jstat <= maxseismo && jstat >= 0) {
                    picks = (NewPick *) realloc(picks, sizeof (NewPick)*(npicks + 1));
                    memset(&picks[npicks], 0, sizeof (NewPick));
                    strcpy(picks[npicks].stname, channel_to_name[jstat]);
                    strcpy(picks[npicks].phase, mprec->p[j].pick.m_phase.phase);
                    strncpy(picks[npicks].quali, mprec->p[j].pick.m_phase.quali, 2);
                    picks[npicks].minute = mprec->p[j].pick.m_phase.minute;
                    picks[npicks].time = mprec->p[j].pick.m_phase.time;
                    picks[npicks].dmy = mprec->p[j].pick.m_phase.dmy;
                    if (picks[npicks].dmy <= 0 || picks[npicks].dmy > maxseismo) {
                        picks[npicks].dmy = jstat + 1;
                    }
                    picks[npicks].weight = mprec->p[j].pick.m_phase.weight;
                    picks[npicks].amplitude = mprec->p[j].pick.m_phase.amplitude;
                    picks[npicks].periode = mprec->p[j].pick.m_phase.period;
                    picks[npicks].amptime = mprec->p[j].pick.m_phase.amptime;
                    strncpy(picks[npicks].ampchan, mprec->p[j].pick.m_phase.ampchan, 4);
                    if (strlen(picks[npicks].ampchan) < 3) {
                        picks[npicks].ampchan[1] = 'H';
                        picks[npicks].ampchan[2] = mprec->p[j].pick.m_phase.ampchan[1];
                    }
                    picks[npicks].ampuse = mprec->p[j].pick.m_phase.ampuse;
                    strncpy(picks[npicks].amptype, mprec->p[j].pick.m_phase.dmy2, 2);
                    if (picks[npicks].amptype[0] != 'W' && picks[npicks].amptype[0] != 'B') picks[npicks].amptype[0] = 0;
                    if (picks[npicks].quali[0] == 'I') {
                        picks[npicks].tmin = -tmin;
                        picks[npicks].tmax = tmin;
                    }
                    if (picks[npicks].quali[0] == 'E') {
                        picks[npicks].tmin = -tmin * 5;
                        picks[npicks].tmax = tmin * 5;
                    }
                    if (picks[npicks].quali[0] == 'Q' || picks[npicks].quali[0] == ' ') {
                        picks[npicks].tmin = -tmin * 10;
                        picks[npicks].tmax = tmin * 10;
                    }

                    npicks++;
                }
                ks++;
            }
        }
        nrec = mprec->p_next;
    }
    if ((pick_type & (short) 1) == 0) {
        old2new_picks(picks, npicks);
    }
    *npicksp = picks;
    return (npicks);
}

void DecodeKPaPicks(m_head, file_data, file_length, lpar)
Master *m_head;
unsigned char *file_data;
int file_length;
NewLocPar *lpar;
{
    int i, j, k;
    if (m_head->a_pick_rec > 0) {
        NautoPicks = ReadOldPicks(m_head->a_pick_rec, file_data, file_length, m_head->pick_type, &autoPicks, &locAutoP);
    } else if (m_head->n_pick_rec > 0) {
        NautoPicks = ReadNewPicks(m_head->n_pick_rec, file_data, file_length, m_head->pick_type, &autoPicks, &locAutoP);
        for (i = 0; i < NautoPicks; i++) {
            for (j = 0; j < nstat; j++) {
                if (strcmp(station[j].name, autoPicks[i].stname) == 0) {
                    for (k = 0; k < station[j].nsig; k++) {
                        if (station[j].shead[k].s_comp[2] == 'Z' &&
                                station[j].shead[k].s_comp[1] == autoPicks[i].ampchan[1]) {
                            autoPicks[i].dmy = station[j].trace_nr[k] + 1;
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }

}

void DecodeKPmPicks(m_head, file_data, file_length)
Master *m_head;
unsigned char *file_data;
int file_length;
{
    if (m_head->m_pick_rec > 0) {
        NmanualPicks = ReadNewPicks(m_head->m_pick_rec, file_data, file_length, m_head->pick_type, &manualPicks, &locManualP);
    }
}

void DecodeKPamLoc(m_head, file_data, file_length)
Master *m_head;
unsigned char *file_data;
int file_length;
{
    int baddr;
    autoLoc.locsum[0] = 0;
    manuLoc.locsum[0] = 0;
    if (m_head->loc_rec > 0) {
        baddr = (m_head->loc_rec - 1) * 512;
        memcpy(&location, &file_data[baddr], sizeof (kpLocRec));
        if (location.auto_loc[0] != 0 && location.auto_loc[0] != ' ') {
            strncpy(autoLoc.locsum, location.auto_loc, 80);
            strncpy(autoLoc.region, location.auto_reg, 80);
            autoLoc.locsum[80] = 0;
            autoLoc.region[80] = 0;
        }
        if (location.manu_loc[0] != 0 && location.manu_loc[0] != ' ') {
            strncpy(manuLoc.locsum, location.manu_loc, 80);
            strncpy(manuLoc.region, location.manu_reg, 80);
            manuLoc.locsum[80] = 0;
            manuLoc.region[80] = 0;
        }
        if (strncmp((char *) & location.fix_mb, "    ", 4) != 0) {
            mag_fixed.fix_ml = location.fix_ml;
            mag_fixed.fix_mb = location.fix_mb;
            mag_fixed.fix_ms = location.fix_ms;
            strncpy(mag_fixed.fix_ms_str, location.fix_ms_str, 4);
            mag_fixed.fix_ms_str[4] = 0;
        }
    }
}

void DecodeKPFilter(m_head, file_data, file_length)
Master *m_head;
unsigned char *file_data;
int file_length;
{
    int baddr;
    int i, n;
    NewFilter FilterPar;
    if (m_head->filter_rec > 0) {
        baddr = (m_head->filter_rec - 1) * 512;
        memcpy(&frec, &file_data[baddr], sizeof (filt_rec));
        n = 0;
        i = 0;
        while (frec.manu_f[i].typ[0] != ' ') {
            manufil = (NewFilter *) realloc(manufil, sizeof (NewFilter)*(n + 1));
            memset(&manufil[n], 0, sizeof (NewFilter));
            strcpy(manufil[n].typ, frec.manu_f[i].typ);
            manufil[n].zerophase = frec.manu_f[i].zerophase;
            manufil[n].low = frec.manu_f[i].low;
            manufil[n].high = frec.manu_f[i].high;
            manufil[n].order = frec.manu_f[i].order;
            i++;
            n++;
        }
        Nmanufil = n;
        n = 0;
        i = 0;
        while (frec.auto_f[i].f_typ[0] != ' ') {
            if (frec.auto_f[i].f_typ[0] == 'T') {
                strcpy(FilterPar.typ, "LP");
                FilterPar.high = 1.0 / frec.auto_f[i].f_corner;
                n++;
            };
            if (frec.auto_f[i].f_typ[0] == 'H') {
                strcpy(FilterPar.typ, "HP");
                FilterPar.low = 1.0 / frec.auto_f[i].f_corner;
                n++;
            }
            FilterPar.zerophase = False;
            FilterPar.order = 2;
            i++;
        }
        if (n >= 2) {
            strcpy(FilterPar.typ, "BP");
            n /= 2;
        }
        if (n != 0) {
            autofil = (NewFilter *) realloc(autofil, sizeof (NewFilter) * n);
            autofil[0] = FilterPar;
            Nautofil = 1;
        }
    }
}

void DecodeKPsignalHead(m_head, mr, nmr, file_data, file_length)
Master *m_head;
R256 *mr;
int nmr;
unsigned char *file_data;
int file_length;
{
    kpSignal *shead;
    char *wave, name[8], comp[8];
    int i, j, l, k, ks, nm, ns, n1;
    int baddr, comp_bit, comp_mask;

    time_minimum = m_head->starttime;
    time_maximum = m_head->endtime;
    sec_minimum = m_head->startsec;
    sec_maximum = m_head->endsec;

    ks = m_head->whereseismo - 1;
    nm = 0;
    s_freq_max = 0;
    ns = m_head->nseismo;
    maxseismo = ns;
    for (i = 0; i < ns; i++) {
        if (ks == 255) {
            ks = 0;
            nm++;
        }
        baddr = (mr[nm][ks] - 1) * 512; /*byte address of signal header */
        shead = (kpSignal *) malloc(sizeof (kpSignal));
        memcpy(shead, &file_data[baddr], sizeof (kpSignal));
        build_name_chanid(shead->s_name, name, comp);
        comp_bit = AddComp2Table(comp);
        comp_mask = 1 << comp_bit;
        n1 = nstat;
        for (j = 0; j < nstat; j++) {
            if (strcmp(name, station[j].name) == 0 && strlen(name) == strlen(station[j].name)) {
                n1 = j;
                break;
            }
        }
        memset(&station[n1].shead[station[n1].nsig], 0, sizeof (SignalHead));
        strcpy(station[n1].shead[station[n1].nsig].s_name, name);
        strcpy(station[n1].shead[station[n1].nsig].s_comp, comp);
        j = SignalHeadFromKP(shead, n1, station[n1].nsig);
        TryRecordingSystem(station[n1].shead[station[n1].nsig].s_resol,
                station[n1].shead[station[n1].nsig].s_volts,
                station[n1].shead[station[n1].nsig].s_recsys,
                &station[n1].shead[station[n1].nsig].s_saturation);
        wave = (char *) malloc(shead->s_nblocks * 512);
        memcpy(wave, &file_data[baddr + 512], shead->s_nblocks * 512);
        j = station[n1].nsig;
        station[n1].wave[j] = wave;
        station[n1].available[j] = comp_mask;
        station[n1].comp_tab[j] = comp_mask;
        strcpy(station[n1].chanid[j], comp);
        if (n1 == nstat) {
            strcpy(station[n1].name, name);
            nstat++;
        }
        if (s_freq_max < station[n1].shead[j].s_freq) s_freq_max = station[n1].shead[j].s_freq;
        if (station[n1].shead[j].s_resol > 0) {
            l = 1;
            for (k = 0; k < station[n1].shead[j].s_resol - 1; k++) l *= 2;
            station[n1].shead[station[n1].nsig].s_saturation = (float) l;
            if (station[n1].shead[j].s_resol < 17) station[n1].shead[j].s_saturation *= 0.75;
        }

        station[n1].trace_nr[j] = i;
        station[n1].nsig++;
        channel_to_station[i] = n1;
        channel_to_compnr[i] = j;
        strcpy(channel_to_name[i], name);
        ks++;
    }
}

void TryRecordingSystem(resol, volts, recsys, sat)
int resol;
float volts;
char *recsys;
float *sat;
{
    int i;

    for (i = 0; i < NrecSys; i++) {
        if (resol == recSys[i].resol) {
            if (fabs(fabs(volts) - recSys[i].volts) < 0.001) {
                strcpy(recsys, recSys[i].recsys);
                *sat = recSys[i].clip;
                break;
            }
        }
    }
}

int SignalHeadFromKP(sh, ns, nsig)
kpSignal *sh;
int ns, nsig;
{
    station[ns].shead[nsig].s_lat = sh->s_lat;
    station[ns].shead[nsig].s_lon = sh->s_lon;
    station[ns].shead[nsig].s_elev = sh->s_elev;
    station[ns].shead[nsig].s_burrial = sh->s_burrial;
    station[ns].shead[nsig].s_start_min = sh->s_start_min;
    station[ns].shead[nsig].s_start_sec = sh->s_start_sec;
    station[ns].shead[nsig].s_end_min = sh->s_end_min;
    station[ns].shead[nsig].s_end_sec = sh->s_end_sec;
    station[ns].shead[nsig].s_freq = sh->s_freq;
    station[ns].shead[nsig].s_nsampl = sh->s_nsampl;
    station[ns].shead[nsig].s_resol = (int) sh->s_resol;
    station[ns].shead[nsig].s_position = (int) sh->s_position;
    station[ns].shead[nsig].s_nbytes = (int) sh->s_nbytes;
    station[ns].shead[nsig].s_nblocks = (int) sh->s_nblocks;
    station[ns].shead[nsig].s_status = (int) sh->s_status;
    station[ns].shead[nsig].s_ndiff = (int) sh->s_ndiff;
    station[ns].shead[nsig].s_n_compress = sh->s_n_compress;
    station[ns].shead[nsig].s_time_base = (int) sh->s_time_base;
    station[ns].shead[nsig].s_gps_sync = (int) sh->s_gps_sync;
    station[ns].shead[nsig].s_volts = sh->s_volts;
    station[ns].shead[nsig].s_gain = sh->s_gain;
    station[ns].shead[nsig].s_damp = sh->s_damp;
    station[ns].shead[nsig].s_const = sh->s_const;
    station[ns].shead[nsig].s_period = sh->s_period;
    station[ns].shead[nsig].s_checksm = sh->s_checksm;
    station[ns].shead[nsig].s_delay = sh->s_delay;
    station[ns].shead[nsig].s_equip = (int) sh->s_equip;
    station[ns].shead[nsig].s_z_sync = (int) sh->s_z_sync;
    station[ns].shead[nsig].s_z_typ = (int) sh->s_z_typ;
    strncpy(station[ns].shead[nsig].s_compress, sh->s_compress, 4);
    station[ns].shead[nsig].s_compress[4] = 0;
    if (strcmp(station[ns].shead[nsig].s_compress, "CMP6") == 0) {
        strcpy(station[ns].shead[nsig].s_compress, "CM6");
    }
    strncpy(station[ns].shead[nsig].s_typ, sh->s_typ, 6);
    station[ns].shead[nsig].s_typ[6] = 0;
    if (sh->s_nsampl > MaxSample) MaxSample = sh->s_nsampl;
    return (0);
}

int SaveKPFilter(file, filter, what)
char *file;
NewFilter *filter;
char *what;
{
    FILE *fp_data;
    int i, k, l, nbytes;
    Master m_head;

    if (strlen(what) == 0) return (1);
    if ((fp_data = fopen(file, "rb+")) == NULL) {
        l = fprintf(stderr, "unable to open file: %s\n", file);
        return (1);
    }
    /****************************************
     * read master record:
     ****************************************/
    nbytes = fread(&m_head, sizeof (Master), 1, fp_data);
    if (NrEvent > 10) {
        i = NrEvent / 10000;
        k = NrEvent - i * 10000;
        m_head.ev_nr = (unsigned short) k;
    }
    if (m_head.filter_rec != 0) {
        i = fseek(fp_data, (long) (m_head.filter_rec - 1)*512, SEEK_SET);
        nbytes = fread(&frec, sizeof (filt_rec), 1, fp_data);
    } else {
        m_head.filter_rec = m_head.lastrec + 1;
        m_head.lastrec++;
        strcpy(frec.manu_f[0].typ, "  ");
    }
    if (strcmp(what, "auto") == 0) {
        strcpy(frec.auto_f[0].f_typ, "HP2");
        frec.auto_f[0].f_corner = 1.0 / filter->low;
        strcpy(frec.auto_f[1].f_typ, "TP2");
        frec.auto_f[1].f_corner = 1.0 / filter->high;
    }
    i = fseek(fp_data, (long) (m_head.filter_rec - 1)*512, SEEK_SET);
    i = fwrite(&frec, sizeof (filt_rec), 1, fp_data);
    i = fseek(fp_data, 0L, SEEK_SET);
    i = fwrite(&m_head, sizeof (Master), 1, fp_data);
    fclose(fp_data);
    return (0);
}

int WritePicksToKPFile(sfile, what, pick, npick, locp)
char *sfile;
char *what;
NewPick *pick;
int npick;
NewLocPar locp;
{
    int rec, mpr_in, mpr_out, i, j, k, l, m, is, nchan;
    unsigned short pr, recnr[200];
    long file_position;
    char chan[4];
    FILE *fp_data = NULL;
    size_t nbytes;
    mpick_rec *pmpr = NULL;
    Master m_head;

    if ((fp_data = fopen(sfile, "rb+")) == NULL) {
        l = fprintf(stderr, "open error for file %s\n", sfile);
        return (1);
    }
    /* read master record */
    nbytes = fread(&m_head, sizeof (Master), 1, fp_data);
    if (NrEvent > 10) {
        i = NrEvent / 10000;
        k = NrEvent - i * 10000;
        m_head.ev_nr = (unsigned short) k;
    }
    mpr_in = 0;
    mpr_out = 0;
    pr = 0;
    for (i = 0; i < strlen(what); i++) what[i] = toupper(what[i]);
    if (npick <= 0) {
        if (what[0] == 'T') m_head.nmx_tr_rec = 0;
        if (what[0] == 'A') m_head.n_pick_rec = 0;
        if (what[0] == 'M') m_head.m_pick_rec = 0;
    } else {

        if (what[0] == 'T') {
            pr = m_head.nmx_tr_rec;
        }
        if (what[0] == 'A') {
            pr = m_head.n_pick_rec;
            m_head.a_pick_rec = 0;
        }
        if (what[0] == 'M') {
            pr = m_head.m_pick_rec;
            m_head.pick_type = 3;
        }
        /* read existing triggers */
        memset(recnr, 0, 200 * sizeof (short));
        if (pr != 0) {
            rec = pr;
            pmpr = (mpick_rec *) malloc(sizeof (mpick_rec));
            memset(pmpr, 0, sizeof (mpick_rec));
            while (rec != 0) {
                file_position = (long) (rec - 1) * 512;
                fseek(fp_data, file_position, SEEK_SET);
                pmpr = (mpick_rec *) realloc(pmpr, sizeof (mpick_rec)*(mpr_in + 1));
                memset(&pmpr[mpr_in], 0, sizeof (mpick_rec));
                nbytes = fread(&pmpr[mpr_in], sizeof (mpick_rec), 1, fp_data);
                recnr[mpr_in] = rec;
                rec = pmpr[mpr_in].p_next;
                mpr_in++;
            }
        }
        mpr_out = (npick + 10) / 10;

        if (mpr_out > mpr_in) {
            pmpr = (mpick_rec *) realloc(pmpr, sizeof (mpick_rec)*(mpr_out));
        }
        memset(pmpr, 0, sizeof (mpick_rec)*(mpr_out));

        k = sortPicks(pick, npick);

        k = 0;
        l = 0;

        pmpr[k].p[l].par.n_picks = npick;
        pmpr[k].p[l].par.st_list = toascii(locp.st_list + 48);
        pmpr[k].p[l].par.model = toascii(locp.model + 48);
        pmpr[k].p[l].par.lat = locp.lat;
        pmpr[k].p[l].par.lon = locp.lon;
        pmpr[k].p[l].par.depth = locp.depth;
        pmpr[k].p[l].par.xnear = locp.xnear;
        pmpr[k].p[l].par.xfar = locp.xfar;
        pmpr[k].p[l].par.pos = locp.pos;
        pmpr[k].p[l].par.dray = locp.dray;
        pmpr[k].p[l].par.daz = locp.daz;
        strncpy(pmpr[k].p[l].par.type, locp.type, 4);
        pmpr[k].p[l].par.fix = toascii(locp.fix + 48);
        pmpr[k].p[l].par.swt = toascii(locp.swt + 48);
        pmpr[k].p[l].par.AmpFilt = toascii(locp.AmpFilt);
        pmpr[k].p[l].par.NoMag = toascii(locp.NoMag);
        pmpr[k].p[l].par.dumy2 = locp.dumy2;

        l++;
        for (i = 0; i < npick; i++) {
            if (l > 9) {
                l = 0;
                k++;
            }
            is = sort_idx[i];
            pmpr[k].p[l].pick.chan = 0;
            pmpr[k].p[l].pick.dummy = 0;
            strcpy(chan, pick[is].ampchan);
            chan[2] = 'Z';
            /* get channel number */
            for (j = 0; j < nstat; j++) {
                if (strcmp(pick[is].stname, station[j].name) != 0) continue;
                nchan = 0;
                for (m = 0; m < station[j].nsig; m++) {
                    if (strcmp(chan, station[j].chanid[m]) == 0) {
                        nchan = station[j].trace_nr[m] + 1;
                        break;
                    }
                }
                if (nchan == 0) {
                    for (m = 0; m < station[j].nsig; m++) {
                        if (chan[2] == station[j].chanid[m][2]) {
                            if (chan[0] == 0) strcpy(pick[is].ampchan, station[j].chanid[m]);
                            nchan = station[j].trace_nr[m] + 1;
                            break;
                        }
                    }
                }
                pmpr[k].p[l].pick.chan = nchan;
                pmpr[k].p[l].pick.dummy = 0;
                break;
            }
            pmpr[k].p[l].pick.m_phase.minute = pick[is].minute;
            pmpr[k].p[l].pick.m_phase.time = pick[is].time;
            strncpy(pmpr[k].p[l].pick.m_phase.phase, pick[is].phase, 8);
            strncpy(pmpr[k].p[l].pick.m_phase.quali, pick[is].quali, 2);
            pmpr[k].p[l].pick.m_phase.weight = pick[is].weight;
            pmpr[k].p[l].pick.m_phase.dmy = pick[is].dmy;
            pmpr[k].p[l].pick.m_phase.amplitude = pick[is].amplitude;
            strncpy(pmpr[k].p[l].pick.m_phase.ampchan, pick[is].ampchan, 4);
            pmpr[k].p[l].pick.m_phase.period = pick[is].periode;
            pmpr[k].p[l].pick.m_phase.amptime = pick[is].amptime;
            pmpr[k].p[l].pick.m_phase.ampuse = pick[is].ampuse;
            strncpy(pmpr[k].p[l].pick.m_phase.dmy2, pick[is].amptype, 2);
            l++;
        }
        /* determine record numbers */
        if (pr == 0) { /* add all new records */
            pr = m_head.lastrec;
            for (i = 0; i < mpr_out; i++) {
                pr++;
                recnr[i] = pr;
                recnr[i + 1] = 0;
            }
            m_head.lastrec = pr;
        } else { /* reuse old records */
            if (mpr_out > mpr_in) {
                pr = m_head.lastrec;
                for (i = mpr_in; i < mpr_out; i++) {
                    pr++;
                    recnr[i] = pr;
                    recnr[i + 1] = 0;
                }
                m_head.lastrec = pr;
            }
            if (mpr_out < mpr_in) {
                recnr[mpr_out] = 0;
            }
        }
        if (what[0] == 'T') {
            m_head.nmx_tr_rec = recnr[0];
        }
        if (what[0] == 'A') {
            m_head.n_pick_rec = recnr[0];
            m_head.a_pick_rec = 0;
        }
        if (what[0] == 'M') {
            m_head.m_pick_rec = recnr[0];
            m_head.pick_type = 3;
        }
        /* now store triggers at specifiyed record position */
        for (i = 0; i < mpr_out; i++) {
            pmpr[i].p_next = recnr[i + 1];
            file_position = (long) (recnr[i] - 1) * 512;
            fseek(fp_data, file_position, SEEK_SET);
            nbytes = fwrite(&pmpr[i], sizeof (mpick_rec), 1, fp_data);
        }
    }
    /* re-write master record */
    /* now re-write master record */
    rewind(fp_data);
    nbytes = fwrite(&m_head, sizeof (Master), 1, fp_data);
    fclose(fp_data);
    return (0);
}

void cparm(argc, argv, param, npar, kget)
int argc;
char *argv[];
OptString *param;
int npar;
int *kget;
{
    int i, k, l;

    for (i = 0; i < nopt; i++) opt[i].set[0] = 0;
    i = 1;
    *kget = 0;
    while (i < argc) {
        for (k = 0; k < nopt; k++) {
            if (strcmp(argv[i], opt[k].def) == 0) {
                if (opt[k].def[0] == '+') {
                    i++;
                    strcpy(opt[k].set, argv[i]);
                } else
                    strcpy(opt[k].set, opt[k].def);
                goto next;
            }
        }
        if (argv[i][0] == '-' || argv[i][0] == '+') {
            l = printf("unknown option %s \n", argv[i]);
            i++;
            continue;
        }
        if (*kget >= npar) {
            l = printf("too many arguments: %s, ignored", argv[i]);
        } else {
            strcpy(param[*kget], argv[i]);
            *kget = *kget + 1;
        };
next:
        i++;
    }

}

int ccheck_opt(which)
char *which;
{
    int i;
    int ret = 1;

    for (i = 0; i < nopt; i++) {
        if (strcmp(opt[i].def, which) == 0) {
            if (opt[i].set[0] != 0) ret = 0;
            break;
        }
    }
    return (ret);
}

char *cget_opt(which)
char *which;
{
    int i;
    char *ret = 0;

    for (i = 0; i < nopt; i++) {
        if (strcmp(opt[i].def, which) == 0) {
            if (opt[i].set[0] != 0) {
                ret = &(opt[i].set[0]);
            };
            break;
        }
    }
    return (ret);
}

void cset_opt(which, string)
char *which;
char *string;
{
    int i;

    for (i = 0; i < nopt; i++) {
        if (strcmp(opt[i].def, which) == 0) {
            strcpy(opt[i].set, string);
            break;
        }
    }
}

/* c- version of juliam and juliam4 */
int juliam(yr, mo, dy, hr, mn)
short *yr, *mo, *dy, *hr, *mn;
{
    int yy, mm, dd, hh, nn;
    yy = *yr;
    mm = *mo;
    dd = *dy;
    hh = *hr;
    nn = *mn;
    return ( juliam4(&yy, &mm, &dd, &hh, &nn));
}

int juliam4(yr, mo, dy, hr, mn)
int *yr, *mo, *dy, *hr, *mn;
{
    int kmo[] = {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    //int leap=1;
    int ky, km, kd, kh, kn, ky4, ky1, ky0, kl, l, minute;
    ky = *yr;
    km = *mo;
    kd = *dy;
    kh = *hr;
    kn = *mn;
    if (km < 0) km = 0;
    minute = ky * 365;
    kd += kmo[km];
    ky4 = ky / 4;
    ky1 = ky / 100;
    ky0 = ky / 1000;
    kl = ky4 - ky1 + ky0;
    if ((ky4 * 4 == ky) && (ky1 * 100 != ky || ky0 * 1000 == ky)) l = 1;
    else l = 0;
    if (l != 0 && km < 3) kl--;
    if (ky == 0 && km == 0) kl = 0;
    minute += kd + kl;
    minute *= 24;
    minute += kh;
    minute *= 60;
    minute += kn;
    return (minute);
}

void datum(julmin, yr, mo, dy, hr, mn)
int *julmin;
short *yr, *mo, *dy, *hr, *mn;
{
    int iyr, imo, idy, ihr, imn;
    datum4(julmin, &iyr, &imo, &idy, &ihr, &imn);
    *yr = iyr;
    *mo = imo;
    *dy = idy;
    *hr = ihr;
    *mn = imn;
}

void datum4(julmin, yr, mo, dy, hr, mn)
int *julmin, *yr, *mo, *dy, *hr, *mn;
{
    int kmo[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int jyr, jmo, jhr, jmn;
    int i, k, kh, id, l, ld, jyr4, jyrh, jyrt;

    k = *julmin / 60;
    jmn = *julmin - k * 60;
    kh = k / 24;
    jhr = k - kh * 24;
    jyr = kh / 365;
    id = -1;
    while (id <= 0) {
        id = kh - jyr * 365;
        l = 0;
        jyr4 = jyr / 4;
        jyrh = jyr / 100;
        jyrt = jyr / 1000;
        ld = jyr4 - jyrh + jyrt;
        if (jyr4 * 4 == jyr && (jyrh * 100 != jyr || jyrt * 1000 == jyr)) l = 1;
        id = id - ld + l;
        jyr--;
    }
    jyr++;
    kmo[2] = 28 + l;
    for (i = 1; i < 13; i++) {
        id -= kmo[i];
        jmo = i;
        if (id <= 0) break;
    }
    *yr = jyr;
    *mo = jmo;
    *dy = id + kmo[jmo];
    *hr = jhr;
    *mn = jmn;
}

/*============================================
  additional routines
 =============================================*/
void build_name_chanid(s0, name, comp)
char *s0, *name, *comp;
{
    int l1, l2;
    char r1[9], s2[9], name8[9];
    strncpy(name8, s0, 8);
    name8[8] = 0;
    csplitstring(name8, r1, s2, ' ');
    l1 = strlen(r1);
    l2 = strlen(s2);
    if (l2 < 3) { /* old chan id */
        if (s2[0] == 'S') { /* short period signals */
            comp[0] = 'S';
            comp[1] = 'H';
            switch (s2[1]) {
                case 'S':
                { /* low gain channels */
                    comp[1] = 'L';
                    comp[2] = 'Z';
                    break;
                }
                case 'M':
                { /* medium gain channels */
                    comp[1] = 'M';
                    comp[2] = 'Z';
                    break;
                }
                case 'G':
                { /* gain-ranging channels */
                    comp[1] = 'R';
                    comp[2] = 'Z';
                    break;
                }
                default:
                { /* high gain channels */
                    comp[1] = 'H';
                    comp[2] = s2[1];
                }
            }
        }
        if (s2[0] == 'B') { /* broad band data */
            comp[0] = 'B';
            comp[1] = 'H';
            comp[2] = s2[1];
        }
        if (s2[0] == 'H') { /* highfrequency broad band data */
            comp[0] = 'H';
            comp[1] = 'H';
            comp[2] = s2[1];
        }
        if (s2[0] == 'L') { /* long period data */
            comp[0] = 'L';
            comp[1] = 'H';
            comp[2] = s2[1];
        }
        if (s2[0] == 'A') { /* accelerometers */
            comp[0] = 'S';
            comp[1] = 'G';
            comp[2] = s2[1];
        }

        /* Changed by hrm 13.12.95
        if (!strncmp(s2,"TC",2))
        {
           if (!strncmp(r1,"INT",3))
           {
              strcpy(comp,"AI ");
           }
           else
           {
              strcpy(comp,"AD ");
           }
        }    end hrm */

        /* Changed by m.b. 25.2.97 */
        comp[3] = 0;
        if (!strncmp(s2, "TC", 2)) {
            comp[0] = 'A';
            comp[1] = r1[0];
            comp[2] = ' ';
            comp[3] = 0;
        } /*  end m.b. */

        if (l2 == 0) { /* this is a timecode or 5-char station name */
            comp[0] = 'A';
            if (l1 > 4) {
                if (strncmp(&r1[4], "TC", 2) == 0) {
                    r1[4] = 0;
                    strcpy(comp, "AD ");
                } else {
                    strncpy(comp, &r1[5], 3);
                    r1[5] = 0;
                    if (strlen(comp) < 3) {
                        comp[2] = comp[1];
                        if (comp[2] != 'S') comp[1] = 'H';
                        else {
                            comp[1] = 'L';
                            comp[2] = 'Z';
                        }
                    }
                }
            } else {
                comp[0] = 'A';
                if (strncmp(r1, "DCF", 3) == 0) comp[1] = 'D';
                if (strncmp(r1, "INT", 3) == 0) comp[1] = 'I';
                comp[2] = ' ';
            }
        }
        comp[3] = 0;
        strcpy(name, r1);
    } else {
        strcpy(name, r1);
        strcpy(comp, s2);
    }
}

int csplitstring(s1, result, s2, c)
char *s1, *result, *s2, c;
/******************************************************* 
   split a string 's1' into result (first part) and
                            's2' (remaining part)
   where the terminating character for is 'c'
   
example: strcpy(s1,"beispiel a-b;c");
         csplitstring("beispiel a-b;c",result,s2,' ');
      --> result: "beispiel"
          s2:     "a-b;c"
         csplitstring("beispiel a-b;c",result,s2,'-');
      --> result: "beispiel a"
          s2:     "b;c"
 *******************************************************/

{
    int i, j, l, k, m;

    l = strlen(s1);
    m = 0;
    for (i = 0; i < l; i++) {
        if (s1[i] == c) {
            result[i] = 0;
            s2[0] = 0;
            k = i + 1;
            for (j = k; j < l; j++) {
                if (s1[j] == c) k++;
                else break;
            }
            strcpy(s2, &s1[k]);
            return (0);
        } else result[i] = s1[i];
        m++;
    }
    s2[0] = 0;
    result[m] = 0;
    return (0);
}

int remove_leading(string, c)
char *string, c;
{
    int i, l, k, n;
    l = strlen(string);
    i = 0;
    while (string[i] == c) i++;
    for (k = i, n = 0; k < l; k++, n++) string[n] = string[k];
    string[n] = 0;
   return(0);
}

char* nfs_mount(package, filename)
char *package;
char *filename;
{
    FILE *sedf;
    char LocalHost[20], RunningHost[20];
    int i, k;

    /* Opening the Package file                                            */
    strcpy(filename, "/signals/info/");
    strcat(filename, package);

    sedf = fopen(filename, "r");
    if (sedf == NULL) {
        strcpy(filename, "/nfs/");
        strcat(filename, package);
        strcat(filename, "/");
        return (filename);
    }
    i = fscanf(sedf, "%s", RunningHost);

    k = gethostname(LocalHost, 20);
    if (strcmp(LocalHost, RunningHost) == 0) {
        strcpy(filename, "/");
    } else {
        if (strlen(RunningHost) > 0) {
            strcpy(filename, "/nfs/");
            strcat(filename, RunningHost);
            strcat(filename, "/");
        } else {
            filename[0] = 0;
        }
    }
    fclose(sedf);
    return (filename);
}

/* 
 * trimlencc expands a string to a new lenth by appending blanks.
 * Caution !!!!!!!!!
 * Be sure, that enough memory is allocated for the new string length
 * trimlencc dont perform any checks!!!
 */
void trimlencc(i_string, max_len, bl_app)
char *i_string; /* in- and output string */
int max_len;
int bl_app;
{
    int a;
    int b;
    int c;

    b = strlen(i_string);
    if (b < max_len) {
        memset(&i_string[b], 32, (size_t) (max_len - b));
        c = b + 1;
    } else c = max_len;
    i_string[max_len] = 0;
    a = 0;
    c--;
    for (a = 0; a < c; a++) {
        if (i_string[a] != 32) b = a;
    };

    /* Place the new terminator */
    i_string[b + bl_app + 1] = 0;
}

void trimlen0c(s1, s2)
char *s1;
char *s2;
{
    int i, l;
    l = strlen(s1);
    strcpy(s2, s1);
    i = 0;
    while (i < l) {
        if (s2[i] == ' ') {
            s2[i] = 0;
            break;
        }
        i++;
    }
}

void trimlenc(s1)
char *s1;
{
    int l;
    l = strlen(s1) - 1;
    while (l >= 0) {
        if (s1[l] != ' ') break;
        s1[l] = 0;
        l--;
    }

}

void trim(s1, s2)
char *s1, *s2;
{
    int i, l, k1, k2;
    l = strlen(s1) - 1;
    i = 0;
    k1 = 0;
    k2 = l;
    while (i < l) {
        if (s1[i] != ' ') {
            k1 = i;
            break;
        }
        i++;
    }
    while (k2 >= 0) {
        if (s1[k2] != ' ') {
            break;
        }
        k2--;
    }
    i = k2 - k1 + 1;
    if (i > 0) {
        strncpy(s2, &s1[k1], i);
    } else {
        i = 0;
    }
    s2[i] = 0;
}

void expandString(s, len, c)
char *s, c;
int len;
{
    int i;
    i = strlen(s);
    while (i < len) {
        s[i] = c;
        i++;
    }
    s[len] = 0;
}

FSize *FindFiles(directory, mask, fn, nn)
char *directory;
char *mask;
FSize *fn;
int *nn;

{
    /* use the routine scandir(3c) to access and sort files alphabetically
       see example inHP-UX Refenece Volume 2
     */

    int num_entries, i, k, l, nm;
    //int year, month, day, hour, min;
    struct dirent **namelist;

    nm = strlen(mask);
    k = *nn;
    if ((num_entries = scandir(directory, &namelist, NULL, alphasort)) < 0) {
        l = fprintf(stderr, "unexpected error\n");
        return (fn);
    }
    if (num_entries) {
        for (i = 0; i < num_entries; i++) {
            if (strncmp(namelist[i]->d_name, mask, nm) == 0) {
                fn = (FSize *) realloc(fn, sizeof (FSize)*(k + 1));
                memset(&fn[k], 0, sizeof (FSize));
                l = sprintf(fn[k], "%s/%s", directory, namelist[i]->d_name);
                k++;
            }
            free(namelist[i]);
        }
        free(namelist);
    }
    *nn = k;
    return (fn);
}

void LimitFiles(fullpath, nl)
char *fullpath;
int nl;
{
    FSize *flist = NULL;
    char clean[256], ext[256], dir[256], bna[256], fn[256];
    int nrlist, *sidx, ne, i, k, l;
    int *nidx;

    if (nl > 0) {
        strcpy(fn, fullpath);
        strcpy(dir, dirname(fn));
        strcpy(fn, fullpath);
        strcpy(bna, basename(fn));
        nrlist = 0;
        flist = FindFiles(dir, bna, flist, &nrlist);
        if (nrlist > 0) {
            nidx = (int *) malloc(nrlist * sizeof (int));
            sidx = (int *) malloc(nrlist * sizeof (int));
            ne = 0;
            for (i = 0; i < nrlist; i++) {
                LastExt(flist[i], clean, ext);
                k = sscanf(ext, "%d", &nidx[ne]);
                if (k < 1) nidx[ne] = -1;
                sidx[ne] = ne;
                ne++;
            }
            if (ne > 1) {
                bubble_sortI(nidx, sidx, ne);
                if (ne <= nl) {
                    l = sprintf(clean, "%s.%d", fullpath, nidx[ne - 1] + 1);
                    l = rename(flist[sidx[ne - 1]], clean);
                }
                for (i = ne - 2; i > 0; i--) {
                    l = rename(flist[sidx[i]], flist[sidx[i + 1]]);
                }
            }
            strcpy(clean, fullpath);
            strcat(clean, ".1");
            l = rename(fullpath, clean);
            free(flist);
            free(nidx);
            free(sidx);
        }
    }
}

FILE *OpenLimiter(file, ext, mode, lim)
char *file;
char *ext;
char *mode;
int lim;
{
    char ofile[256];
    strcpy(ofile, file);
    if (strlen(ext) > 0) {
        strcat(ofile, ".");
        strcat(ofile, ext);
    }
    if (mode[0] == 'w') {
        LimitFiles(ofile, lim);
    }
    return (fopen(ofile, mode));
}

void GetStationList(file)
char *file;
{
    FILE *fp = NULL; /* Pointer to open file   */
    char Line[1024];
    int l;

    gse_nst = 0;
    if (strlen(file) == 0) {
        strcpy(file, "/usr/local/lib/sed_stations.GSE_SED");
    }
    if ((fp = fopen(file, "r")) != NULL) {
        fgets(Line, 1024, fp);
        while (!feof(fp)) {
            if (Line[0] == '(') {
                fgets(Line, 1024, fp);
                continue;
            }
            l = strlen(Line);
            gse_stn = (GseStatList *) realloc(gse_stn, sizeof (GseStatList)*(gse_nst + 1));
            memset(&gse_stn[gse_nst], 0, sizeof (GseStatList));
            sscanf(Line, "%s %s %f %f %f",
                    gse_stn[gse_nst].name,
                    gse_stn[gse_nst].comp,
                    &gse_stn[gse_nst].lat,
                    &gse_stn[gse_nst].lon,
                    &gse_stn[gse_nst].elev);
            if (Line[40] != ' ' && l > 40) strncpy(gse_stn[gse_nst].ondate, &Line[40], 10);
            if (Line[51] != ' ' && l > 51) strncpy(gse_stn[gse_nst].ofdate, &Line[51], 10);
            if (Line[80] != ' ' && l > 80) strncpy(gse_stn[gse_nst].nation, &Line[80], 3);
            if (l > 119) sscanf(&Line[119], "%d %f", &gse_stn[gse_nst].v_model, &gse_stn[gse_nst].MOHO_depth);
            if (Line[132] != ' ' && l > 132) gse_stn[gse_nst].bulletin = Line[132];
            if (l > 150)
                sscanf(&Line[134], "%f %f %f",
                    &gse_stn[gse_nst].coeff[0],
                    &gse_stn[gse_nst].coeff[1],
                    &gse_stn[gse_nst].coeff[2]);
            gse_nst++;
            fgets(Line, 1024, fp);
        }
        fclose(fp);
    }
}

void Basename(fullpath, path, name)
char *fullpath, *path, *name;
{
    strcpy(str1, fullpath);
    strcpy(path, dirname(str1));
    strcpy(str1, fullpath);
    strcpy(name, basename(str1));
}
