/***********************************************************************
** I give ALATURKA in an "as is" basis. I do not say or imply that    **
** it will be useful for whatever you want to do with it. You can     **
** use ALATURKA free of charge in academic research and teaching.     **
** For any commercial use, you must get a written permission from me. **
**                                                                    **
** (C) Ali Dasdan, 1991-                                              **
***********************************************************************/
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <malloc.h>
# include "../common/defs.h"
# include "defs_sa.h"
# include "../hcommon/random.h"
# include "../hcommon/fileio.h"
# include "../hcommon/readinput.h"
# include "../hcommon/partition.h"
# include "lib.h"
# include "lib_sa.h"

#define DEBUG1 
#undef DEBUG2 

/* Simulated Annealing for Multiple-way Hypergraph Partitioning - based on */
/* Peterson-etal's version. */

/* definitions */
cells_type            cells [MAX_CELLS];
nets_type             nets [MAX_NETS];
corn_type             cnets  [MAX_PINS];
corn_type             ncells [MAX_PINS];
ind_type              pop [MAX_POP];             /* population */
cells_info_type       cells_info [MAX_CELLS];   
                                         /* additional information for cells */
selected_cell_type    scell [1];     /* selected cell type */
selected_cell_type    prev_scell [1];     /* previously selected cell type */
eval_type             *eval;

/* partitioning variables */
int nocells;           /* number of cells */
int noparts;           /* number of partitions */
int nonets;            /* number of nets */
int nopins;            /* number of pins */
int totcellsize;       /* total cell weight of the partition */
int totnetsize;        /* total net weight of the partition */
int cutsize;           /* cutsize of the partition */
int max_gain;          /* max gain of a cell */
int max_cdeg;          /* max density of a cell */
int max_ndeg;          /* max density of a net */
int max_cweight;       /* max cell weight */
int max_nweight;       /* max net weight */
int mov_count;         /* count of total moves */
float off_ratio;       /* alpha in initial partition alg. */

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
int region1;           /* flags to identify the regions during entire run of SA */
int region2;
int region3;

main (argc, argv)
int argc;
char *argv [];
{
  FILE *fp;
  char fname [STR_SIZE];
  long seed;
  int  i, pass_no;

  if (argc < 3) {
    printf ("\nUsage: InputFileName noparts [seed]\n");
    exit (1);
  }  /* if */
  sprintf (fname, "../chip/%s", argv [1]);
  noparts = atoi (argv [2]);                         /* (D) */
  if (argc > 3)   seed = (long) atoi (argv [3]);
  else   seed = (long) -1;
  seed = randomize ((long)  seed);
  printf ("SEED = %d fname = %s\n", seed, fname);
  /* read input hypergraph */
  read_hgraph (fname, &nocells, &nonets, noparts, &nopins,
              &totcellsize, &totnetsize, &max_cdeg, &max_ndeg,
              &max_cweight, &max_nweight,
              cells, nets, cnets, ncells);
  off_ratio = (float) 0.1;   /* used in balance criterion */
  /* create an initial partition */
  create_partition (nocells, noparts, totcellsize, max_cweight, &off_ratio,
                    cells, nets, cnets, &pop [0]);
  printf ("off=%f\n", off_ratio);
  printf ("Initial : Part_no min_size curr_size max_size\n");
  for (i = 0; i < noparts; i++) {
    printf ("II %d %d %d %d\n", i, pop [0].parts [i].pmin_size,
            pop [0].parts [i].pcurr_size, pop [0].parts [i].pmax_size);
  }
  max_gain = max_cdeg * max_nweight;
  /* find initial cutsize */
  cutsize = find_cut_size (nonets, noparts, totnetsize, nets, &pop [0]);
  printf ("Totalsize = %d Initial cutsize = %d\n", totnetsize, cutsize);
  /* compute costs (gains) of cells in hypergraph */
  compute_gains (nocells, noparts, cells, nets, cnets, cells_info, pop [0].chrom);
  /* initialize variables of SA algorithm */
  temperature = 10.0;
  templength = nocells * (noparts - 1); 
  /* the whole run has three different regions */
  region1 = True;    
  region2 = region3 = False;
  pass_no = 0;
  frozen = False;
  while (! frozen)  {
    costsum = 0.0;
    cost2sum = 0.0;
    nouphills = 0;
    nochanges = 0;    
    changed = False;
    for (notrials = 0; notrials < templength; notrials++) {
      /* randomly select a cell to move to a randomly selected part */
      selected = select_cell (nocells, noparts, scell, pop [0].chrom, 
                              cells, pop [0].parts, cells_info);
      if (! selected) {
        printf ("Cannot find a move to select.\n");
        exit (1);
      }   /* if */
      delta = - scell [0].mov_gain;
      if (((float) delta <= 0.0) ||
          (((float) delta > 0.0) && 
           (rand01 () <= (float) exp ((double) -delta / (double) temperature)))) {
#ifdef DEBUG2
      printf ("cell=%d from=%d to=%d gain=%d\n",
              scell [0].mov_cell_no, scell [0].from_part,
              scell [0].to_part, scell [0].mov_gain);
#endif
        /* move the selected cell */
        move_cell (scell, pop [0].chrom, cells, pop [0].parts);
        /* update the costs of the neighbor cells */
        update_gains (scell, cells, nets, cnets, ncells, 
                      cells_info, pop [0].chrom);
        cutsize += delta;
        /* update variables due to the change */
        changed = True;
        nochanges++;
        costsum += (float) cutsize;
        cost2sum += (float) (cutsize * cutsize);
        if ((float) delta > 0.0) nouphills++;
      }   /* if the selected move is accepted */
    }   /* for */
#ifdef DEBUG1
      printf ("changes=%d trials=%d \% accept=%f temperature=%f\n", 
              nochanges, notrials, 
              (100.0 * (float) nochanges / notrials), temperature);
#endif
    /* update variables of SA algorithm */
    pass_no++;
    printf ("pass_no = %d Region=%d%d%d Final cutsize  = %d Check cutsize = %d\n",
             pass_no, region1, region2, region3, cutsize, 
             find_cut_size (nonets, noparts, totnetsize, nets, &pop [0]));
    if (changed) {
      if (region1) {   /* when in region1 of the run */
        costvariance = cost2sum / nochanges - 
                       ((costsum * costsum) / (nochanges * nochanges));
        if (temperature <= 0.0) {
          printf (" Temperature is non-positive.\n");
          exit (1);
        }   /* if */
        costvariance = costvariance / temperature;
        if (costvariance >= 0.05) tempfactor = 1.25;   /* = 1.0 / 0.8 */
        else {
          region1 = False; region2 = True;
        }   /* else */
      }   /* if region1 */
      else if (region2) {   /* when in region2 of the run */
        percentchanges = (float) nochanges / templength;
        if (percentchanges > 0.5) tempfactor = 0.95;
        else {
          region2 = False; region3 = True;
        }   /* if */
      }   /* else if region2 */
      else if (region3) {   /* when in region3 of the run */
        if (nouphills > 0) {
          tempfactor = 0.95;
          templength = nocells * (noparts - 1) * (noparts - 1);
        }   /* if */
        else region3 = False;
      }   /* else if region3 */
      else
        frozen = True;
      temperature = tempfactor * temperature;
    }   /* if changed */
    else {
      printf ("No change during this temperature length.\n");
      exit (1);
    }   /* else if not changed */
  }   /* while not frozen */ 
  /* final output */
  printf ("Final : Part_no min_size curr_size max_size\n");
  for (i = 0; i < noparts; i++) {
    printf ("FF %d %d %d %d\n", i, pop [0].parts [i].pmin_size,
            pop [0].parts [i].pcurr_size, pop [0].parts [i].pmax_size);
  }
  exit (0);
}   /* main */
