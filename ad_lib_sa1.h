#ifndef AD_LIB_PFM_INCLUDED
#define AD_LIB_PFM_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* randomly select a move among (nocells * (noparts - 1)) move directions */
int select_cell(int nocells, 
                int noparts,
                selected_cell_t scell[],
                allele tchrom[],
                cells_t cells[],
                parts_t parts[],
                cells_info_t cells_info[]);

/* find a to part to which the selected cell is to be moved */
/* return -1 if not found */
int find_to_part(int noparts, 
                 int cell_no,
                 int from,
                 cells_t cells[],
                 parts_t parts[]);

/* return true if the selected move is feasible, else false */
int feasible_move(int cell_no,
                  int from,
                  int to,
                  cells_t cells[],
                  parts_t parts[]);

/* move the selected cell */
void move_cell(selected_cell_t scell[],
               allele tchrom[],
               cells_t cells[],
               parts_t parts[]);

/* update gains after a move */
void update_gains(selected_cell_t scell[],
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  corn_t ncells[],
                  cells_info_t cells_info[],
                  allele tchrom[]);

/* update gain of a cell */
void update1(int flag,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             nets_t nets[],
             corn_t ncells[],
             cells_info_t cells_info[]);

/* update gain of a cell */
void update2(int flag,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             nets_t nets[],
             corn_t ncells[],
             cells_info_t cells_info[],
             allele tchrom[]);

/* copy pop structures */
void copy_pop(int nocells,
              int noparts,
              ind_t pop1[], 
              ind_t pop2[]);

/* copy nets structures */
void copy_nets(int nonets,
               int noparts,
               nets_t from_nets[],
               nets_t to_nets[]);

#endif
