#ifndef AD_PRINT_INCLUDED
#define AD_PRINT_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

void print_hgraph(int nocells,
                  int nonets,
                  cells_t cells[],
                  nets_t nets[],
                  corn_t ncells[]);

void print_vars(int nocells,
                int nonets,
                int noparts,
                int nopins,
                int totcellsize,
                int totnetsize,
                int max_cdeg,
                int max_ndeg,
                int max_cweight,
                int max_nweight);

void print_cnets(int nopins, corn_t cnets[]);

void print_ncells(int nopins, corn_t ncells[]);

void print_parts(int nocells, 
                 int noparts,
                 ind_t *ind,
                 allele tchrom[]);

void print_parts_info(int nocells,
                      int noparts,
                      allele chrom[],
                      parts_info_t *parts_info);

void print_cells(int nocells,
                 cells_t cells[],
                 corn_t cnets[]);

void print_nets(int nonets,
                int noparts,
                nets_t nets[],
                corn_t ncells[]);

void print_chrom(int nocells, allele chrom[]);

void print_cells_info(int nocells,
                      int noparts,
                      cells_info_t cells_info[]);

void print_inx(int noparts, partb_t **partb);

void print_buckets(int noparts,
                   int bucketsize,
                   partb_t **partb);

#endif
