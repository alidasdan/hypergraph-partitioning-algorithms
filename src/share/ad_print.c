
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include "ad_defs.h"
#include "ad_print.h"

void print_hgraph(int nocells,
                  int nonets,
                  cells_t cells[],
                  nets_t nets[],
                  corn_t ncells[])
{
    printf("\nHYPERGRAPH\n");
    printf("%d\n%d\n", nocells, nonets);
    for (int i = 0; i < nonets; i++) {
        printf("%d %d ", nets[i].nweight, nets[i].nno_cells);
        for (int j = 0; j < nets[i].nno_cells; j++) 
            printf("%d ", ncells[j + nets[i].celllist].corn_no);
        printf("\n"); 
    }   /* for i */
    for (int i = 0; i < nocells; i++) 
        printf("%d\n", cells[i].cweight); 
}   /* print_hgraph */

void print_vars(int nocells,
                int nonets,
                int noparts,
                int nopins,
                int totcellsize,
                int totnetsize,
                int max_cdeg,
                int max_ndeg,
                int max_cweight,
                int max_nweight)
{
    printf("\nHYPERGRAPH VARIABLES\n");
    printf("nocells= %d\n", nocells);
    printf("noparts= %d\n", noparts);
    printf("nonets= %d\n", nonets);
    printf("nopins= %d\n", nopins);
    printf("totcellsize= %d\n", totcellsize);
    printf("totnetsize= %d\n", totnetsize);
    printf("max_cdeg= %d\n", max_cdeg);
    printf("max_ndeg= %d\n", max_ndeg);
    printf("max_cweight= %d\n", max_cweight);
    printf("max_nweight= %d\n", max_nweight);
}   /* print_vars */

void print_cnets(int nopins, corn_t cnets[])
{
    printf("\nCNETS\n");
    for (int i = 0; i < nopins; i++) {
        printf("%d ", cnets[i].corn_no);
        if ((i % 10) == 0) {
            printf("\n");
        }
    }
    printf("\n");
}   /* print_cnets */

void print_ncells(int nopins, corn_t ncells[])
{
    printf("\nNCELLS\n");
    for (int i = 0; i < nopins; i++) {
        printf("%d ", ncells[i].corn_no);
        if ((i % 10) == 0) {
            printf("\n");
        }
    }
    printf("\n");
}   /* print_ncells */

void print_parts(int nocells, 
                 int noparts,
                 ind_t *ind,
                 allele tchrom[])
{
    printf("\nPARTS\n");

    for (int i = 0; i < noparts; i++) {

        printf("part_no=%d: %d <= %d <= %d * #cells=%d\n", i,
               ind->parts[i].pmin_size,
               ind->parts[i].pcurr_size,
               ind->parts[i].pmax_size,
               ind->parts[i].pmax_cells);

        int k = 0;
        for (int j = 0; j < nocells; j++) {
            if (tchrom[j] == i)  {
                printf("%d ", j);
                k++;
                if (k % 20 == 0) {
                    printf("\n");
                }
            }   /* if */
        }   /* for j */

        printf("\n");

    }   /* for i */
}   /* print_parts */

void print_parts_info(int nocells,
                      int noparts,
                      allele chrom[],
                      parts_info_t *parts_info)
{
    printf("\nPARTS_INFO\n");

    for (int i = 0; i < noparts; i++) {

        printf("part_no=%d: %d <= %d <= %d * #cells=%d\n", i,
               parts_info[i].pmin_size,
               parts_info[i].pcurr_size,
               parts_info[i].pmax_size,
               parts_info[i].pmax_cells);

        for (int j = 0; j < nocells; j++) {
            if (chrom[j] == i) {
                printf("%d ", j);
            }
        }   /* for j */

        printf("\n");

    }   /* for i */
}   /* print_parts_info */

void print_cells(int nocells,
                 cells_t cells[],
                 corn_t cnets[])
{
    printf("\nCELLS\n");

    int totsize = 0;
    for (int i = 0; i < nocells; i++) {

        printf("c=%d #n=%d cw=%d nl=%d ns= ", 
               i, cells[i].cno_nets, cells[i].cweight, cells[i].netlist);

        totsize += cells[i].cweight;

        for (int j = 0; j < cells[i].cno_nets; j++)
            printf("%d ", cnets[cells[i].netlist + j].corn_no);

        printf("\n");

    }   /* for */

    printf("totcellsize = %d\n", totsize);
}   /* print_cells */

void print_nets(int nonets,
                int noparts,
                nets_t nets[],
                corn_t ncells[])
{
    printf("\nNETS\n");

    int totsize = 0;
    for (int i = 0; i < nonets; i++) {

        printf("n=%d #c=%d nw=%d cl=%d p(n)= ", i, nets[i].nno_cells,
               nets[i].nweight, nets[i].celllist);

        totsize += nets[i].nweight;

        for (int j = 0; j < noparts; j++) 
            printf("%d:%d ", j, nets[i].npartdeg[j]);

        printf("\ncs= ");

        for (int j = 0; j < nets[i].nno_cells; j++) 
            printf("%d ", ncells[nets[i].celllist + j].corn_no); 

        printf("\n");

    }   /* for i*/
 
    printf("totnetsize= %d\n", totsize);
}   /* print_nets */

void print_chrom(int nocells, allele chrom[])
{
    printf("\nCHROM:\n");
    for (int i = 0; i < nocells; i++) {
        printf("%d ", chrom[i]);
    }
    printf("\n");
}   /* print_chrom */

void print_cells_info(int nocells,
                      int noparts,
                      cells_info_t cells_info[])
{
    printf("\nCELLS_INFO\n");
    for (int i = 0; i < nocells; i++) {
        printf("c=%d mc=%d l=%d mg= ", i, cells_info[i].mcount, cells_info[i].locked);
        for (int j = 0; j < noparts; j++) 
            printf("%d ", cells_info[i].mgain[j]);
        printf("\n");
    }   /* for i */
}   /* print_cells_info */

void print_inx(int noparts, partb_t **partb)
{
    printf("\nINDICES:\n");
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {
            printf("(%d,%d) max=%d min=%d nobucs=%d\n",
                   i, j, partb[i][j].max_inx, partb[i][j].min_inx, partb[i][j].nobuckets);
        }  /* for j */
    }  /* for i */
}   /* print_inx */

void print_buckets(int noparts,
                   int bucketsize,
                   partb_t **partb)
{
    printf("\nPARTB:\n");

    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {
            for (int k = 0; k < bucketsize; k++) {

                bnode_ptr_t next = partb[i][j].bnode_ptr[k];
                if (next != NULL) {
                    printf("Partb[%d][%d].bnode_ptr[%d]=", i, j, k);
                    while (next != NULL) {
                        if (next->lptr != NULL) {
                            printf("%d<", (next->lptr)->cell_no);
                        }
                        printf("%d>", next->cell_no);
                        if (next->rptr != NULL) {
                            printf("%d ", (next->rptr)->cell_no);
                        }
                        next = next->rptr;
                    }   /* while */

                    printf("\n");

                }   /* if */
            }   /* for k */
        }   /* for j */
    }   /* for i */
}   /* print_buckets */

/* EOF */
