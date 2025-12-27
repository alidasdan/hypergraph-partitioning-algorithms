
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdlib.h>
#include <string.h>
#include <malloc/malloc.h>
#include <errno.h>
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_defs.h"
#include "ad_partition.h"

/* create initial partition */
int create_partition(int nocells,
                     int noparts,
                     int totcellsize,
                     int max_cweight,
                     float *off_ratio,
                     cells_t cells[],
                     nets_t nets[],
                     corn_t cnets[],
                     ind_t *ind)
{
    /* init current size */
    for (int i = 0; i < noparts; i++) {
        ind->parts[i].pratio = 1.0 / noparts;
        ind->parts[i].pcurr_size = 0;   
        ind->parts[i].pmax_cells = 0;
    }   /* for i */

    /* allocate temporary memory */
    int *tparts = (int *) calloc(noparts, sizeof(int));
    if (tparts == NULL) {
        printf("Error: Unable to allocate memory for tparts.\n");
        exit(1);
    }   /* if */

    /* insert cells */
    for (int i = 0; i < nocells; i++) {

        /* find minimum */
        int min_inx = 0;
        int min_size = ind->parts[0].pcurr_size;
        for (int j = 1; j < noparts; j++) {
            if (min_size > ind->parts[j].pcurr_size) {
                min_inx = j;
                min_size = ind->parts[j].pcurr_size;
            }   /* if */
        }   /* for j */

        /* find a new minimum among the same minimums */
        int tcount = -1; /* tparts is used to randomize partitioning */
        for (int j = 0; j < noparts; j++) {
            if (ind->parts[j].pcurr_size == min_size) {
                tcount++; 
                tparts[tcount] = j;
            }   /* if */
        }   /* for j */
        min_inx = tparts[irandom(0, tcount)];

        /* assign cell i to part[min_inx] */
        ind->chrom[i] = min_inx;
        ind->parts[min_inx].pcurr_size += cells[i].cweight; 
        ind->parts[min_inx].pmax_cells++; 

        /* find net info */
        int cnets_inx = cells[i].netlist;
        for (int j = 0; j < cells[i].cno_nets; j++) {
            int net_no = cnets[cnets_inx + j].corn_no;
            nets[net_no].npartdeg[min_inx]++;
        }   /* for j */

    }   /* for i */

    free(tparts);

    /* determine min & max part sizes */
    /* also adjust the min & max part sizes */
    int part_fit = True;
    int max_size = -1;
    while (((*off_ratio) < 1.05) && (part_fit)) {

        for (int i = 0; i < noparts; i++) {
            float part_size = ((float) totcellsize) * ind->parts[i].pratio;
            if (part_size < max_cweight) {
                printf("\nError: Too small part size.\n");
                exit(1);
            }   /* if */
            ind->parts[i].pmax_size = (int) (part_size * (1.0 + (*off_ratio)) + 1.0);
            ind->parts[i].pmin_size = (int) (part_size * (1.0 - (*off_ratio)) - 1.0);
            if (ind->parts[i].pmax_size > max_size) {
                max_size = ind->parts[i].pmax_size;
            }
        }   /* for i */ 

        int i = 0;
        while ((i < noparts) && (part_fit)) {
            if (ind->parts[i].pmax_size < ind->parts[i].pcurr_size)
                part_fit = False;
            i++;
        }   /* while */

        if (part_fit) {
            part_fit = False;
        } else {
            part_fit = True;
            (*off_ratio) += 0.05;
            if ((*off_ratio) > 1.0) {
                printf("\nError: Cannot adjust part sizes.\n");
                exit(1);
            }   /* if */
        }   /* if not fit */

    }   /* while */

    return max_size;
}   /* create_partition */

/* copy partition properties to parts_info for temporary use */
void copy_partition(int noparts,
                    parts_info_t *parts_info,
                    ind_t *ind)
{
    for (int i = 0; i < noparts; i++) {
        parts_info[i].pmax_cells = ind->parts[i].pmax_cells;
        parts_info[i].pmin_size = ind->parts[i].pmin_size;
        parts_info[i].pcurr_size = ind->parts[i].pcurr_size;
        parts_info[i].pmax_size = ind->parts[i].pmax_size;
    }   /* for i */
}   /* copy_partition */

/* read a partition prepared beforehand */
int read_partition(FILE *fp,
                   char *filename,
                   int noparts,
                   cells_t cells[],
                   nets_t nets[],
                   corn_t cnets[],
                   ind_t *ind)
{
    int max_size = -1; 
    open_file(&fp, filename, "r");

    for (int i = 0; i < noparts; i++) {

        int part_no;
        if (fscanf(fp, "%d", &part_no) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", filename, errno, strerror(errno));
        }

        if (fscanf(fp, "%d%d%d%d",
                   &(ind->parts[part_no].pmax_cells),
                   &(ind->parts[part_no].pmin_size),
                   &(ind->parts[part_no].pcurr_size),
                   &(ind->parts[part_no].pmax_size)) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", filename, errno, strerror(errno));        
        }

        if (ind->parts[part_no].pmax_size > max_size)
            max_size = ind->parts[part_no].pmax_size;

        for (int j = 0; j < ind->parts[part_no].pmax_cells; j++) {

            int cell_no;
            if (fscanf(fp, "%d", &cell_no) == EOF) {
                printf("Error: Cannot read from %s: errno= %d error= %s\n", filename, errno, strerror(errno));
            }
            ind->chrom[cell_no] = part_no; 

            /* find net info */
            int cnets_inx = cells[cell_no].netlist;
            for (int k = 0; k < cells[cell_no].cno_nets; k++) {
                int net_no = cnets[cnets_inx + k].corn_no;
                nets[net_no].npartdeg[part_no]++;
            }   /* for k */

        }   /* for j */

    }   /* for i */

    close_file(&fp);

    return max_size;
}   /* read_partition */

/* write a partition */
void write_partition(FILE *fp,
                     char *filename,
                     int nocells,
                     int noparts,
                     ind_t *ind)
{
    open_file(&fp, filename, "w");

    for (int i = 0; i < noparts; i++) {
        fprintf(fp, "%d %d %d %d %d\n", i,
                ind->parts[i].pmax_cells,
                ind->parts[i].pmin_size,
                ind->parts[i].pcurr_size,
                ind->parts[i].pmax_size);
        for (int j = 0; j < nocells; j++)
            if (ind->chrom[j] == i)
                fprintf(fp, "%d ", j);
        fprintf(fp, "\n");
    }   /* for i */

    close_file(&fp);
}   /* read_partition */

/* EOF */
