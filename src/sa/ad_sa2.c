
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
   based on the paper: Peterson et al., Neural networks and
   NP-complete problems; a performance study of the graph bisectioning
   problem, Complex Systems 2, 1989. */

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
    float costsum;         /* sum of cutsizes during a temperature length */
    float cost2sum;        /* sum of squares of cutsizes during a temperature length */
    float percentchanges;  /* percentage of accepted moves (= changes) */
    float costvariance;    /* variance of cutsizes (= cost) */ 
    float temperature;     /* temperature in SA alg. */
    float tempfactor;      /* cooling ratio */
    int templength;        /* number of runs until equilibrium at a certain temperature */
    int nochanges;         /* number of accepted states */
    int notrials;          /* number of tried states */
    int selected;          /* set if a feasible move is found */
    int changed;           /* set if a state is accepted */ 
    int delta;             /* cost difference between the last two states */
    int nouphills;         /* number of uphill moves */
    int frozen;            /* set if the system is frozen - after all three regions */
    int region2;           /* flags to identify the regions during entire run of SA */
    int region3;
    int region4;
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

#if DEBUG
    printf("Totalsize = %d Initial cutsize = %d\n", totnetsize, cutsize);
#endif

    /* compute costs (gains) of cells in hypergraph */
    compute_gains(nocells, noparts, cells, nets, cnets, cells_info, pop[0].chrom);

    /* initialize the current champion solution */
    best_cutsize = cutsize;
    copy_pop(nocells, noparts, pop, best_pop);
    copy_nets(nonets, noparts, nets, best_nets);

    /* initialize variables of SA algorithm */
    temperature = 10.0;
    templength = nocells * (noparts - 1); 

    /* the whole run has three different regions: region1 is the initialization. */
    region2 = True;    
    region3 = region4 = False;
    pass_no = 0;

    /* while not yet frozen */
    frozen = False;
    while (! frozen)  {

        costsum = 0.0;
        cost2sum = 0.0;
        nouphills = 0;
        nochanges = 0;    
        changed = False;

        /* perform the following exploration templength times */
        for (notrials = 0; notrials < templength; notrials++) {

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
                printf ("cell=%d from=%d to=%d gain=%d\n",
                        scell[0].mov_cell_no, scell[0].from_part,
                        scell[0].to_part, scell[0].mov_gain);
#endif

                /* move the selected cell */
                move_cell(scell, pop[0].chrom, cells, pop[0].parts);

                /* update the costs of the neighbor cells */
                update_gains(scell, cells, nets, cnets, ncells, 
                             cells_info, pop[0].chrom);
                cutsize += delta;

                /* update the current champion solution */
                if (cutsize < best_cutsize) {
                    best_cutsize = cutsize;
                    copy_pop(nocells, noparts, pop, best_pop);
                    copy_nets(nonets, noparts, nets, best_nets);
                }

                /* update variables due to the change */
                changed = True;
                nochanges++;
                costsum += (float) cutsize;
                cost2sum += (float) (cutsize * cutsize);
                if ((float) delta > 0.0) {
                    nouphills++;
                }

            }   /* if the selected move is accepted */
        }   /* for */

#ifdef DEBUG1
        printf("changes=%d trials=%d \% accept=%f temperature=%f\n", 
               nochanges, notrials, 
               (100.0 * (float) nochanges / notrials), temperature);
#endif

        /* update variables of SA algorithm */
        pass_no++;

#ifdef DEBUG1
        printf("pass_no = %d Region=%d%d%d Final cutsize  = %d Check cutsize = %d\n",
               pass_no, region2, region3, region4, cutsize, 
               find_cut_size(nonets, noparts, totnetsize, best_nets, &best_pop[0]));
#endif

        if (changed) {
            tempfactor = 0.95;

            if (region2) {   /* when in region2 of the run: heating up */
                costvariance = cost2sum / nochanges - 
                    ((costsum * costsum) / (nochanges * nochanges));
                costvariance = costvariance / temperature;
                if (costvariance >= 0.05) {
                    tempfactor = 1.25;   /* = 1.0 / 0.8 */
                } else {
                    region2 = False; region3 = True;
                }
            } else if (region3) {   /* when in region3 of the run: cooling */
                percentchanges = (float) nochanges / templength;
                if (percentchanges > 0.5) { 
                    tempfactor = 0.95;
                } else {
                    region3 = False; region4 = True;
                } 
            } else if (region4) {   /* when in region4 of the run: slow cooling */
                if (nouphills > 0) {
                    tempfactor = 0.95;
                    templength = nocells * (noparts - 1) * (noparts - 1);
                } else {
                    region4 = False;
                }
            } else {
                frozen = True;
            }

            /* reduce the temperature */
            temperature = tempfactor * temperature;
        } else {
            printf("Warning: No change during this temperature length.\n");
            exit(1);
        }   /* else if not changed */
    }   /* while not frozen */ 

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
}   /* main-sa2 */

/* EOF */
