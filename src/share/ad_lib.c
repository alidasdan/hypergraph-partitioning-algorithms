
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include <stdlib.h>
#include "ad_defs.h"
#include "ad_lib.h"

/* initialize all bucket indices and pointers */
void init_buckets(int noparts, 
                  int bucketsize,
                  partb_t partb[][noparts - 1])
{
    /* init partition bucket indices */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {

            partb[i][j].max_inx = partb[i][j].min_inx = -1;
            partb[i][j].nobuckets = 0;

            /* init partb bucket pointers */
            for (int k = 0; k < bucketsize; k++) {
                partb[i][j].bnode_ptr[k] = NULL;
            }

        }   /* for j */
    }   /* for i */
}   /* init_buckets */
 
/* map part no such that home_part is excluded */
int map_part_no(int dest_part, int home_part)
{
    if (dest_part < home_part) {
        return (dest_part);
    } else if (dest_part > home_part) {
        return (dest_part - 1);
    } else {
        printf("Error: Unexpected inputs\n");
        exit(1);
    }
}   /* map_part_no */
 
/* compute move gain from home_part to dest_part */
int calculate_gain(int cell_no, 
                   int home_part, 
                   int dest_part,
                   cells_info_t cells_info[])
{
    int mov_gain = cells_info[cell_no].mgain[dest_part] - 
        cells_info[cell_no].mgain[home_part];
    return mov_gain;
}   /* calculate_gain */

/* compute gains of all cells and place them into cells_info */
void compute_gains(int nocells, 
                   int noparts,
                   cells_t cells[],
                   nets_t nets[],
                   corn_t cnets[],
                   cells_info_t cells_info[],
                   allele tchrom[])
{
    for (int cell_no = 0; cell_no < nocells; cell_no++) {

        /* initialize cells_info */
        cells_info[cell_no].locked = False;
        cells_info[cell_no].mcount = 0;

        /* find info about the cell */
        int net_ptr = cells[cell_no].netlist;
        cells[cell_no].cno_inets = 0; 

        /* initialize external & internal costs */
        for (int j = 0; j < noparts; j++) {
            cells_info[cell_no].mgain[j] = 0;
        }
        int part_no = tchrom[cell_no];

        /* for each net on the cell cell_no */
        for (int j = 0; j < cells[cell_no].cno_nets; j++) {

            int net_no = cnets[net_ptr].corn_no;
            int net_weight = nets[net_no].nweight;

            /* if only one cell in part_no */
            /* else if all cells in part_no */
            if (nets[net_no].npartdeg[part_no] == 1) {

                /* find the part in which all other cells of net lie */
                int found = False;
                int tpart_no = 0;
                while ((tpart_no < noparts) && (!found)) {
                    if ((tpart_no != part_no) && 
                        (nets[net_no].npartdeg[tpart_no] == (nets[net_no].nno_cells - 1))) {
                        found = True;
                        cells_info[cell_no].mgain[tpart_no] += net_weight;
                    }   /* if tpart_no */
                    tpart_no++;
                }   /* while */

            } else if (nets[net_no].npartdeg[part_no] == nets[net_no].nno_cells) {
                /* update internal cost */
                cells_info[cell_no].mgain[part_no] += net_weight;
            }   /* else */
            net_ptr++;

        }   /* for all nets j on cell cell_no */

    }   /* for all cells cell_no */
}   /* compute_gains */

/* compute gains of all cells and place them into cells_info */
void compute_gains2(int nocells, 
                    int noparts,
                    cells_t cells[],
                    nets_t nets[],
                    corn_t cnets[],
                    cells_info_t cells_info[],
                    allele tchrom[],
                    nets_info_t nets_info[])
{
    for (int cell_no = 0; cell_no < nocells; cell_no++) {

        /* initialize cells_info */
        cells_info[cell_no].locked = False;
        cells_info[cell_no].mcount = 0;

        /* find info about the cell */
        cells[cell_no].cno_inets = 0; 
        int net_ptr = cells[cell_no].netlist;

        /* initialize external & internal costs */
        for (int j = 0; j < noparts; j++) {
            cells_info[cell_no].mgain[j] = 0;
        }
        int part_no = tchrom[cell_no];

        /* for each net on cell_no */
        for (int j = 0; j < cells[cell_no].cno_nets; j++) {

            int net_no = cnets[net_ptr].corn_no;
            int net_weight = nets[net_no].nweight;

            /* if only one cell in part_no */
            /* else if all cells in part_no */
            if (nets_info[net_no].npartdeg[part_no] == 1) {

                /* find the part in which all other cells of net lie */
                int found = False;
                int tpart_no = 0;
                while ((tpart_no < noparts) && (!found)) {
                    if ((tpart_no != part_no) && 
                        (nets_info[net_no].npartdeg[tpart_no] == 
                         (nets[net_no].nno_cells - 1))) {
                        found = True;
                        cells_info[cell_no].mgain[tpart_no] += net_weight;
                    }   /* if tpart_no */
                    tpart_no++;
                }   /* while */

            } else if (nets_info[net_no].npartdeg[part_no] == nets[net_no].nno_cells) {
                /* update internal cost */
                cells_info[cell_no].mgain[part_no] += net_weight;
            }   /* else */
            net_ptr++;

        }   /* for all nets j on cell cell_no */

    }   /* for all cells cell_no */
}   /* compute_gains2 */

/* free all allocated nodes */
void free_nodes(int noparts, 
                int bucketsize,
                partb_t partb[][noparts - 1])
{
    /* delete nodes connected to partb */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {

            for (int k = 0; k < bucketsize; k++) {
                bnode_ptr_t next = partb[i][j].bnode_ptr[k];
                partb[i][j].bnode_ptr[k] = NULL;
                while (next != NULL) {
                    bnode_ptr_t prev = next;
                    next = next->rptr;
                    free(prev);
                }   /* while */
            }   /* for k */

            partb[i][j].max_inx = -1;
            partb[i][j].min_inx = -1;
            partb[i][j].nobuckets = 0;
        }   /* for j */
    }   /* for i */
}   /* free_nodes */

/* count number of bucket nodes */
void number_nodes(int noparts, 
                  int bucketsize, 
                  int *npartb,
                  partb_t partb[][noparts - 1])
{
    *npartb = 0;

    /* count nodes connected to partb */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {
            for (int k = 0; k < bucketsize; k++) {
                bnode_ptr_t next = partb[i][j].bnode_ptr[k];
                while (next != NULL) {
                    next = next->rptr;
                    (*npartb)++;
                }   /* while */
            }   /* for k */
        }   /* for j */
    }   /* for i */
}   /* number_nodes */

/* find set of cells to be actually moved */
int find_move_set(mcells_t mcells[],
                  int msize,
                  int *max_mcells_inx)
{
    int max_gain_sum = 0;
    *max_mcells_inx = -1;
    int gain_sum = 0;
    for (int i = 0; i < msize; i++) {
        gain_sum += mcells[i].mgain;
        if (gain_sum > max_gain_sum) {
            *max_mcells_inx = i;
            max_gain_sum = gain_sum;
        }   /* if */
    }   /* for i */

    return max_gain_sum;
}   /* find_move_set */

/* move cells actually - permanently */
int move_cells(int wflag,
               int nocells, 
               int msize,
               mcells_t mcells[],
               int max_mcells_inx, 
               int cutsize, 
               int *glob_inx,
               ind_t *ind,
               cells_t cells[],
               nets_t nets[],
               corn_t cnets[])
{
    /* used only under wflag == True */
    int tcutsize, fcutsize;
    tcutsize = fcutsize = cutsize;

    if (wflag == True) {
        if ((*glob_inx) == 0) {
            printf("%d %d %d\n", *glob_inx, tcutsize, fcutsize);
            (*glob_inx)++;
        }   /* if */
    }   /* if */

    int cut_gain = 0;
    for (int i = 0; i <= max_mcells_inx; i++) {

        if (wflag == True) {
            tcutsize -= mcells[i].mgain;
            fcutsize = tcutsize;
            printf("%d %d %d\n", *glob_inx, tcutsize, fcutsize);
            (*glob_inx)++;
        }    /* if wflag */

        ind->chrom[mcells[i].cell_no] = mcells[i].to;
        cut_gain += mcells[i].mgain;

        /* update partition size limits */
        ind->parts[mcells[i].from].pmax_cells--;
        ind->parts[mcells[i].from].pcurr_size -= cells[mcells[i].cell_no].cweight;
        ind->parts[mcells[i].to].pmax_cells++;
        ind->parts[mcells[i].to].pcurr_size += cells[mcells[i].cell_no].cweight;

        /* find the nets on the cell moved and update its npartdeg */
        int cnets_inx = cells[mcells[i].cell_no].netlist;
        for (int j = 0; j < cells[mcells[i].cell_no].cno_nets; j++) {
            int net_no = cnets[cnets_inx + j].corn_no;
            nets[net_no].npartdeg[mcells[i].from]--;
            nets[net_no].npartdeg[mcells[i].to]++;
        }   /* for j */

    }   /* for i */

    for (int i = max_mcells_inx + 1; i < msize; i++) {
        if (wflag == True) {
            tcutsize -= mcells[i].mgain;
            printf("%d %d %d\n", *glob_inx, tcutsize, fcutsize);
            (*glob_inx)++;
        }    /* if wflag */
    }   /* for i */

    for (int i = msize; i < nocells; i++) {
        if (wflag == True) {
            printf("%d %d %d\n", *glob_inx, fcutsize, fcutsize);
            (*glob_inx)++;
        }    /* if wflag */
    }   /* for i */

    return cut_gain;
}   /* move_cells */

/* finds cut size of a given partition - used for control */
int find_cut_size(int nonets, 
                  int noparts, 
                  int totnetsize,
                  nets_t nets[],
                  ind_t *ind)
{
    ind->incost = 0;
    for (int i = 0; i < nonets; i++) {

        int found = False;
        int j = 0;
        while ((j < noparts) && (!found)) {
            if (nets[i].npartdeg[j] == nets[i].nno_cells) {
                found = True;
            }
            j++;
        }   /* while */
        if (found) {
            ind->incost += nets[i].nweight;
        }
    }   /* for i */

    return (totnetsize - ind->incost);
}   /* find_cut_size */

/* save nets info in nets_info */
void copy_nets_info(int nonets, 
                    int noparts,
                    nets_t nets[],
                    nets_info_t nets_info[])
{
    for (int i = 0; i < nonets; i++) {
        for (int j = 0; j < noparts; j++) {
            nets_info[i].npartdeg[j] = nets[i].npartdeg[j];
        }   /* for j */
    }   /* for i */
}   /* copy_nets_info */

/* EOF */
