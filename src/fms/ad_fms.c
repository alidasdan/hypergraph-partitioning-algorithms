
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <assert.h> 
#include <stdio.h>
#include <string.h>
#include <malloc/malloc.h>
#include "ad_defs.h"
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_readinput.h"
#include "ad_partition.h"
#include "ad_print.h"
#include "ad_bucketio.h"
#include "ad_lib.h"
#include "ad_lib_fms.h"
 
/* FOR SANCHIS' VERSION OF MULTI-WAY PARTITIONING */
/* Also mentioned as the SN algorithm */
/* Direct multi-way partitioning.
   Locking is used.
   Cells are moved wrt their gains
*/

int main(int argc, char *argv[])
{
    /* definitions */
    int nocells;           /* number of cells */
    int nonets;            /* number of nets */
    int nopins;            /* number of pins */
    int noparts;           /* number of partitions */
    int totcellsize;       /* total cell weight of the partition */
    int totnetsize;        /* total net weight of the partition */
    int cutsize;           /* cutsize of the partition */
    int max_gain;          /* max gain of a cell */
    int max_cdeg;          /* max density of a cell */
    int max_ndeg;          /* max density of a net */
    int max_cweight;       /* max cell weight */
    int max_nweight;       /* max net weight */
    int bucketsize;        /* max size of a bucket array */
    int msize;             /* index to mcells */

    if (argc < 3) {
        printf("\nUsage: %s InputFileName NoParts [Seed]\n", argv[0]);
        exit(1);
    }  /* if */

    char fname[STR_SIZE];
    sprintf(fname, "%s", argv[1]);

    noparts = atoi(argv[2]);                         

    long seed;
    if (argc > 3) {
        seed = (long) atoi(argv[3]);
    } else {
        seed = -1;
    }
    seed = randomize((long) seed);
    printf("SEED = %ld fname = %s\n", seed, fname);

    read_hgraph_size(fname, &nocells, &nonets, &nopins);

    /* alloc memory for all data structures */
    cells_t *cells = (cells_t *) calloc(nocells, sizeof(cells_t));
    assert(cells != NULL);
    cells_info_t *cells_info = (cells_info_t *) calloc(nocells, sizeof(cells_info_t));
    assert(cells_info != NULL);
    for (int i = 0; i < nocells; i++) {
        cells_info[i].mgain = (int *) calloc(noparts, sizeof(int));
        assert(cells_info[i].mgain != NULL);
        cells_info[i].partb_ptr = (bnode_ptr_t *) calloc(noparts - 1, sizeof(bnode_ptr_t));
        assert(cells_info[i].partb_ptr != NULL);
        cells_info[i].partb_gain_inx = (int *) calloc(noparts - 1, sizeof(int));
        assert(cells_info[i].partb_gain_inx != NULL);
    }

    nets_t *nets = (nets_t *) calloc(nonets, sizeof(nets_t));
    assert(nets != NULL);
    nets_info_t *nets_info = (nets_info_t *) calloc(nonets, sizeof(nets_info_t));
    assert(nets_info != NULL);
    for (int i = 0; i < nonets; i++) {
        nets[i].npartdeg = (int *) calloc(noparts, sizeof(int));
        assert(nets[i].npartdeg != NULL);
        nets_info[i].npartdeg = (int *) calloc(noparts, sizeof(int));
        assert(nets_info[i].npartdeg != NULL);
    }

    /* cells of nets */
    corn_t *cnets = (corn_t *) calloc(nopins, sizeof(corn_t));
    assert(cnets != NULL);
    /* nets of cells */
    corn_t *ncells = (corn_t *) calloc(nopins, sizeof(corn_t));
    assert(ncells != NULL);

    /* partition buckets */
    partb_t partb[noparts][noparts - 1];  
    parts_info_t parts_info[noparts]; 

    /* population (w/ one individual!) */
    ind_t pop[MAX_POP];             
    for (int i = 0; i < MAX_POP; i++) {
        pop[i].chrom = (allele *) calloc(nocells, sizeof(allele));
        assert(pop[i].chrom != NULL);
        pop[i].parts = (parts_t *) calloc(noparts, sizeof(parts_t));
        assert(pop[i].parts != NULL);
    }

    /* selected cell */
    selected_cell_t scell[1];     

    /* moved cells */
    mcells_t *mcells = (mcells_t *) calloc(nocells, sizeof(mcells_t));
    assert(mcells != NULL);
 
    read_hgraph(fname, nocells, nonets, nopins, noparts, 
                &totcellsize, &totnetsize, &max_cdeg, &max_ndeg,
                &max_cweight, &max_nweight,
                cells, nets, cnets, ncells);

    max_gain = max_cdeg * max_nweight;
    bucketsize = 2 * max_gain + 1;

    /* alloc memory (statically if possible) */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            partb[i][j].bnode_ptr = (bnode_ptr_t *) calloc(bucketsize, sizeof(bnode_ptr_t));
        }
    }

    float off_ratio = (float) 0.1;
    create_partition(nocells, noparts, totcellsize, max_cweight, &off_ratio, 
                     cells, nets, cnets, &pop[0]);

#ifdef DEBUG
    printf("off=%f\n", off_ratio);
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0]. parts[i].pmax_size);
    }
#endif

    init_buckets(noparts, bucketsize, partb);
    cutsize = find_cut_size(nonets, noparts, totnetsize, nets, &pop[0]);

#ifdef DEBUG
    printf("Totalsize = %d Initial cutsize = %d\n", totnetsize, cutsize);
#endif

    int gain_sum;
    int no_iter = 0;
    int glob_inx = 0;
    do {

        copy_partition(noparts, parts_info, &pop[0]);

        copy_nets_info(nonets, noparts, nets, nets_info);

        compute_gains(nocells, noparts, cells, nets, cnets, 
                      cells_info, pop[0].chrom);

        create_buckets(nocells, noparts, max_gain, pop[0].chrom,
                       partb, cells_info);

        msize = 0;

        int nlocked = 0;
        do {
            int move_possible = select_cell(noparts, scell, parts_info, cells,
                                            partb, cells_info);

            delete_partb_nodes_of_cell(noparts, scell[0].mov_cell_no,
                                       scell[0].from_part, partb, cells_info);

            /* lock cell */
            cells_info[scell[0].mov_cell_no].locked = True;
            if (move_possible == True) {
                move_cell(mcells, msize, scell);
                msize++;
                update_gains(noparts, max_gain, scell, 
                             cells, nets, cnets, ncells, nets_info, 
                             partb, cells_info, pop[0].chrom);
            }   /* if */
            nlocked++;

        } while (nlocked < nocells);

        int max_mcells_inx;
        gain_sum = find_move_set(mcells, msize, &max_mcells_inx);

#ifdef DEBUG
        printf("gain_sum=%d max_mcells_inx=%d msize = %d\n",
               gain_sum, max_mcells_inx, msize);
#endif

        if (gain_sum > 0) {
            int cut_gain = move_cells(False, nocells, msize, mcells, max_mcells_inx, 
                                      cutsize, &glob_inx, &pop[0], cells, nets, cnets);
            cutsize -= cut_gain;
        }   /* if */
        no_iter++;

#ifdef DEBUG
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", 
               no_iter, cutsize, find_cut_size(nonets, noparts, totnetsize, nets, &pop[0]));
#endif

    } while ((gain_sum > 0) && (cutsize > 0) && (no_iter < NO_ITERATIONS));

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", no_iter, 
           cutsize, find_cut_size(nonets, noparts, totnetsize, nets, &pop[0]));

    free_nodes(noparts, bucketsize, partb);

#ifdef DEBUG
    printf("Final : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("FF %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    /* free memory for all data structures */
    free(cells);
    for (int i = 0; i < nocells; i++) {
        free(cells_info[i].mgain);
        free(cells_info[i].partb_ptr);
        free(cells_info[i].partb_gain_inx);
    }
    free(cells_info);

    for (int i = 0; i < nonets; i++) {
        free(nets[i].npartdeg);
        free(nets_info[i].npartdeg);
    }
    free(nets);
    free(nets_info);

    free(cnets);
    free(ncells);

    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            free(partb[i][j].bnode_ptr);
        }
    }

    for (int i = 0; i < MAX_POP; i++) {
        free(pop[i].chrom);
        free(pop[i].parts);
    }

    free(mcells);

    return (0);
}   /* main-fms */

/* EOF */
