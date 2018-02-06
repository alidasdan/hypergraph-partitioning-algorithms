
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc/malloc.h>
#include "ad_defs.h"
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_readinput.h"
#include "ad_partition.h"
#include "ad_lib.h"
#include "ad_lib_sa.h"

/* Simulated Annealing for Multiple-way Hypergraph Partitioning -
   based on the paper: Johnson et al., Optimization by Simulated
   Annealing, Operations Research, 37(6), 1989. */

int main(int argc, char *argv[])
{
    /* partitioning variables */
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

    /* SA algorithm variables */
    float temperature;     /* temperature in SA alg. */
    float tempfactor;      /* cooling ratio */
    float minpercent;      /* percentage of accepted moves */
    float cutoff;          /* a speedup option - see the paper */
    int delta;             /* cost difference between the last two states */
    int freezelimit;       /* when freezecount reaches this limit, the system is frozen */
    int freezecount;      
    int nochanges;         /* number of accepted states */
    int changeslimit;      /* limit on max number of accepted changes */
    int notrials;          /* number of tried states */
    int trialslimit;       /* limit on max number of tried states */
    int sizefactor;        /* propartionality constant for temperature length */
    int neigh_size;        /* neighborhood size */
    int selected;          /* set if a feasible move is found */
    int changed;           /* set if a change occurs in best_cutsize */
    int same;              /* count the number of times cutsize remains the same */
    int samecount;         /* limit on the counter "same" */ 
    int prev_cutsize;      /* previous cutsize value - used with "same" */
    int pass_no;           /* pass number */

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
        cells_info[i].partb_ptr = NULL;
        cells_info[i].partb_gain_inx = NULL;
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

    /* properties of the best solution */
    int best_cutsize;

    nets_t *best_nets = (nets_t *) calloc(nonets, sizeof(nets_t));
    assert(best_nets != NULL);
    for (int i = 0; i < nonets; i++) {
        best_nets[i].npartdeg = (int *) calloc(noparts, sizeof(int));
        assert(best_nets[i].npartdeg != NULL);
    }

    ind_t best_pop[MAX_POP];
    for (int i = 0; i < MAX_POP; i++) {
        best_pop[i].chrom = (allele *) calloc(nocells, sizeof(allele));
        assert(best_pop[i].chrom != NULL);
        best_pop[i].parts = (parts_t *) calloc(noparts, sizeof(parts_t));
        assert(best_pop[i].parts != NULL);
    }

    read_hgraph(fname, nocells, nonets, nopins, noparts, 
                &totcellsize, &totnetsize, &max_cdeg, &max_ndeg,
                &max_cweight, &max_nweight,
                cells, nets, cnets, ncells);

    /* create the initial partition or solution */
    float off_ratio = (float) 0.1;   /* alpha in initial partitioning */
    create_partition(nocells, noparts, totcellsize, max_cweight, &off_ratio,
                     cells, nets, cnets, &pop[0]);

#ifdef DEBUG1
    printf("off=%f\n", off_ratio);
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    max_gain = max_cdeg * max_nweight;

    /* find initial cutsize */
    cutsize = find_cut_size(nonets, noparts, totnetsize, nets, &pop[0]);
#ifdef DEBUG1
    printf("Totalsize = %d Initial cutsize = %d\n", totnetsize, cutsize);
#endif

    /* compute costs (gains) of cells in hypergraph */
    compute_gains(nocells, noparts, cells, nets, cnets, cells_info, pop[0].chrom);

    /* initialize the current champion solution */
    best_cutsize = cutsize;
    copy_pop(nocells, noparts, pop, best_pop);
    copy_nets(nonets, noparts, nets, best_nets);

    /* initialize variables of SA algorithm */
    changed = False;
    neigh_size = nocells * (noparts - 1);
    sizefactor = 16;      /* 16 in paper */
    minpercent = 0.02;    /* 0.02 in paper */
    cutoff = 1.0;         /* 1.0 in paper */
    temperature = 10.0;   /* 0.6 - not given in paper */
    tempfactor = 0.95;    /* 0.95 in paper */
    freezelimit = 5;      /* 5 in paper */
    freezecount = 0;
    pass_no = 0;
    same = 0;
    samecount = nocells / 2; 
    prev_cutsize = -1;
    trialslimit = sizefactor * neigh_size;
    changeslimit = (int) (cutoff * sizefactor * neigh_size);

    /* while not yet frozen */
    while (freezecount < freezelimit)  {

        /* perform the following exploration trialslimit times */
        nochanges = notrials = 0;
        while ((notrials < trialslimit) && (nochanges < changeslimit)) {
            notrials++;

            /* randomly select a cell to move to a randomly selected part */
            selected = select_cell(nocells, noparts, scell, pop[0].chrom, 
                                   cells, pop[0].parts, cells_info);
            if (! selected) {
                printf("Error: Cannot find a move to select.\n");
                exit(1);
            }   /* if */

            /* delta is the change in the cost function */
            delta = -scell[0].mov_gain;

            /* if delta is negative, this change is a downhill move so
               accept the new solution */
            /* if delta is positive, this change is an uphill move so
               accept with an ever decreasing probability that depends
               on delta and the current temperature */
            if (((float) delta <= 0.0) ||
                (((float) delta > 0.0) && 
                 (rand01() <= (float) exp((double) -delta / (double) temperature)))) {

#ifdef DEBUG2
                printf("cell=%d from=%d to=%d gain=%d\n",
                       scell[0].mov_cell_no, scell[0].from_part,
                       scell[0].to_part, scell[0].mov_gain);
#endif

                nochanges++;

                /* move the selected cell */
                move_cell(scell, pop[0].chrom, cells, pop[0].parts);

                /* update the costs of the neighbor cells */
                update_gains(scell, cells, nets, cnets, ncells, 
                             cells_info, pop[0].chrom);
                cutsize += delta;

                /* update the current champion solution */
                if (cutsize < best_cutsize) {
                    changed = True;
                    best_cutsize = cutsize;
                    copy_pop(nocells, noparts, pop, best_pop);
                    copy_nets(nonets, noparts, nets, best_nets);
                }   /* if */

            }   /* if the selected move is accepted */
        }   /* while */

#ifdef DEBUG1
        printf("changes=%d trials=%d accept=%f temperature=%f\n", 
               nochanges, notrials, 
               (100.0 * (float) nochanges / notrials), temperature);
#endif

        /* reduce the temperature and update the other variables of SA
           algorithm */
        temperature = tempfactor * temperature; 
        if (changed) {
            freezecount = 0;
        }
        if (((float) nochanges / notrials) < minpercent) {
            freezecount++;
        }
        pass_no++;
        changed = False;
        if (cutsize == prev_cutsize) {
            same++;
        } else {
            same = 0;
            prev_cutsize = cutsize;
        }
        /* if has seen the same solution enough number of times, it is
           time to quit */
        if (same >= samecount) { 
            freezecount = freezelimit;  /* exit */
        }

#ifdef DEBUG1
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n",
               pass_no, best_cutsize, 
               find_cut_size(nonets, noparts, totnetsize, 
                             best_nets, &best_pop[0]));
#endif

    }   /* while */ 

#ifdef DEBUG1
    printf("Why : same=%d accept rate=%f\n", same, (float) nochanges / notrials);
#endif

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n",
           pass_no, best_cutsize, 
           find_cut_size(nonets, noparts, totnetsize, best_nets, &best_pop[0]));

#ifdef DEBUG1
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
    }
    free(cells_info);

    for (int i = 0; i < nonets; i++) {
        free(nets[i].npartdeg);
        free(nets_info[i].npartdeg);
        free(best_nets[i].npartdeg);
    }
    free(nets);
    free(nets_info);
    free(best_nets);

    free(cnets);
    free(ncells);

    for (int i = 0; i < MAX_POP; i++) {
        free(pop[i].chrom);
        free(pop[i].parts);
        free(best_pop[i].chrom);
        free(best_pop[i].parts);
    }

    return (0);
}   /* main-sa1 */

/* EOF */
