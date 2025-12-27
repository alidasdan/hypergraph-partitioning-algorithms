#ifndef AD_LIB_PLM_INCLUDED
#define AD_LIB_PLM_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* map a given mov_gain into index to a bucket array */
int map_gain(int mov_gain, int max_gain);
 
/* fill all bucket arrays */
void create_buckets(int nocells,
                    int noparts,
                    int max_gain,
                    allele chrom[],
                    partb_t **partb,
                    cells_info_t cells_info[]);

/* select a cell to move */
int select_cell(int noparts,
                selected_cell_t scell[],
                parts_info_t *parts_info,
                cells_t cells[],
                partb_t **partb,
                cells_info_t cells_info[]);

/* move selected cell, and save the move in a file */
void move_cell(mcells_t mcells[],
               int msize,
               selected_cell_t scell[],
               allele tchrom[]);

/* update gains after a move */
void update_gains(int noparts,
                  int max_gain,
                  selected_cell_t scell[],
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  corn_t ncells[],
                  nets_info_t nets_info[],
                  partb_t **partb,
                  cells_info_t cells_info[],
                  allele tchrom[]);

/* update gain of a cell */
void update1(int flag,
             int noparts,
             int max_gain,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             nets_t nets[],
             corn_t ncells[],
             partb_t **partb,
             cells_info_t cells_info[],
             allele tchrom[]);

/* update gain of a cell */
void update2(int flag,
             int noparts,
             int max_gain,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             nets_t nets[],
             corn_t ncells[],
             partb_t **partb,
             cells_info_t cells_info[],
             allele tchrom[]);

void create_partb_nodes_of_cell(int noparts,
                                int max_gain,
                                int cell_no,
                                int part_no,
                                partb_t **partb,
                                cells_info_t cells_info[]);

#endif

