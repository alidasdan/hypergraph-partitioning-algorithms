
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "ad_defs.h"
#include "ad_fileio.h"
#include "ad_readinput.h"

/* initialize cells array */
void init_cells(int nocells, cells_t cells[])
{
    for (int i = 0; i < nocells; i++) {
        cells[i].cno_nets = 0;
        cells[i].cno_inets = 0;
        cells[i].netlist = NIL;
    }   /* for */
}   /* init_cells */

/* initialize nets array */
void init_nets(int nonets,
               int noparts,
               nets_t nets[])
{
    for (int i = 0; i < nonets; i++) {
        nets[i].nno_cells = 0;
        nets[i].celllist = NIL;
        for (int j = 0; j < noparts; j++) {
            nets[i].npartdeg[j] = 0;
        }
    }   /* for i */
}   /* init_nets */

/* initialize netlist pointers */
void init_netlist(int nonets,
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  corn_t ncells[])
{
    for (int i = 0; i < nonets; i++) {
        for (int j = 0; j < nets[i].nno_cells; j++) {
            int cell_no = ncells[nets[i].celllist + j].corn_no;
            int cnets_inx = cells[cell_no].netlist + cells[cell_no].cno_inets;
            cnets[cnets_inx].corn_no = i;
            cells[cell_no].cno_inets++;
            if (cells[cell_no].cno_inets > cells[cell_no].cno_nets) {
                printf("Error: Inconsistency in cell_%d degrees.\n", j);
                exit(1);
            }   /* if */
        }   /* for j */
    }   /* for i */
}   /* init_netlist */

void read_hgraph_size(char fname[],
                      int  *nocells,
                      int  *nonets,
                      int  *nopins)
{
    FILE *fp;

    open_file(&fp, fname, "r");

    if (fscanf(fp, "%d%d%d", nocells, nonets, nopins) == EOF) {
        printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
        close_file(&fp);
        exit(1);
    } 

    if ((*nocells < 0) || (*nonets < 0) || (*nopins < 0)) {
        printf("Error: Invalid attributes of graph.\n");
        close_file(&fp);
        exit(1);
    }   /* if */

    close_file(&fp);
}   /* read_hgraph_size */

/* read input hypergraph and construct cells, nets, cnets, & ncells arrays */
void read_hgraph(char fname[],
                 int  nocells,
                 int  nonets,
                 int  nopins,
                 int  noparts,
                 int *totcellsize,
                 int *totnetsize,
                 int *max_cdeg,
                 int *max_ndeg,
                 int *max_cweight,
                 int *max_nweight,
                 cells_t cells[],
                 nets_t  nets[],
                 corn_t  cnets[],
                 corn_t  ncells[])
{
    FILE *fp;
    open_file(&fp, fname, "r");

    /* hgraph size is already read so re-read and discard. */
    int ignore1, ignore2, ignore3;
    if (fscanf(fp, "%d%d%d", &ignore1, &ignore2, &ignore3) == EOF) {
        printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
        close_file(&fp);
        exit(1);
    } 

    /* initialize cells & nets arrays */
    init_cells(nocells, cells);
    init_nets(nonets, noparts, nets);

    /* initialize variables */
    *max_cweight = -1;
    *max_nweight = -1;
    *max_cdeg = -1;
    *max_ndeg = -1;
    *totcellsize = 0;
    *totnetsize = 0;

    /* read nets */
    int ncells_inx = 0;
    for (int i = 0; i < nonets; i++) {
        if (fscanf(fp, "%d%d", &nets[i].nweight, &nets[i].nno_cells) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
            close_file(&fp);
            exit(1);
        }
        (*totnetsize) += nets[i].nweight;
        if (nets[i].nweight > (*max_nweight)) {
            *max_nweight = nets[i].nweight;
        }
        if (nets[i].nno_cells > (*max_ndeg)) {  
            *max_ndeg = nets[i].nno_cells;
        }
        nets[i].celllist = ncells_inx;

        for (int j = 0; j < nets[i].nno_cells; j++) {
            int cell_no;
            if (fscanf(fp, "%d", &cell_no) == EOF) {
                printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
                close_file(&fp);
                exit(1);
            }
            ncells[ncells_inx].corn_no = cell_no;
            ncells_inx++;
            cells[cell_no].cno_nets++;
        }   /* for j */
    }   /* for i */

    /* read  cell weights */
    int cnets_inx = 0;
    for (int i = 0; i < nocells; i++) {
        if (fscanf(fp, "%d", &cells[i].cweight) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
            close_file(&fp);
            exit(1);
        }
        (*totcellsize) += cells[i].cweight;
        if (cells[i].cweight > (*max_cweight))
            *max_cweight = cells[i].cweight;
        if (cells[i].cno_nets > (*max_cdeg))
            *max_cdeg = cells[i].cno_nets;
        cells[i].netlist = cnets_inx;
        cnets_inx += cells[i].cno_nets;
    }   /* for i */

    close_file(&fp);

    /* create netlists */
    init_netlist(nonets, cells, nets, cnets, ncells);
}   /* read_hgraph */

/* EOF */
