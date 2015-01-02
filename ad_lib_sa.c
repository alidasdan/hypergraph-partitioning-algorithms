
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ad_defs.h"
#include "ad_random.h"
#include "ad_bucketio.h"
#include "ad_lib.h"
#include "ad_lib_sa.h"

/* randomly select a move among (nocells * (noparts - 1)) move directions */
int select_cell(int nocells, 
                int noparts,
                selected_cell_t scell[],
                allele tchrom[],
                cells_t cells[],
                parts_t parts[],
                cells_info_t cells_info[])
{
    int cell_no = irandom(0, nocells - 1);
    int selected = False;    /* true if a cell is selected */
    int nochecks = 0;    /* number of cells tried */
    while ((! selected) && (nochecks < nocells)) {
        nochecks++;
        int from = tchrom[cell_no]; 
        int to = find_to_part(noparts, cell_no, from, cells, parts);
        if (to < 0) {  /* no directions for this cell are feasible, so try another */
            cell_no = (cell_no + 1) % nocells;
        } else {
            selected = True;
            scell[0].mov_cell_no = cell_no;
            scell[0].from_part = from;
            scell[0].to_part = to;
            scell[0].mov_gain = calculate_gain(scell[0].mov_cell_no,
                                               scell[0].from_part,
                                               scell[0].to_part,
                                               cells_info);
        }   /* else */
    }   /* while */
    if (nochecks < nocells) {
        return (True);   /* a move is found */
    } else {
        return (False);
    }
}   /* select_cell */

/* find a to part to which the selected cell is to be moved */
/* return -1 if not found */
int find_to_part(int noparts, 
                 int cell_no,
                 int from,
                 cells_t cells[],
                 parts_t parts[])
{
    int to = irandom(0, noparts - 1);
    int nochecks = 0;
    int selected = False;
    while ((! selected) && (nochecks < noparts)) {
        nochecks++;
        if ((to == from) || (! feasible_move(cell_no, from, to, cells, parts))) {
            to = (to + 1) % noparts;
        } else {
            selected = True;
        }
    }   /* while */
    if (selected) {
        return (to);
    } else {
        return (-1);
    }
}   /* find_to_part */

/* return true if the selected move is feasible, else false */
int feasible_move(int cell_no,
                  int from,
                  int to,
                  cells_t cells[],
                  parts_t parts[])
{
    int fpcurr_size = parts[from].pcurr_size;
    int fpmin_size = parts[from].pmin_size;
    int tpcurr_size = parts[to].pcurr_size;
    int tpmax_size = parts[to].pmax_size;
    if ((fpcurr_size >= (fpmin_size + cells[cell_no].cweight)) &&
        ((tpcurr_size + cells[cell_no].cweight) <= tpmax_size)) {
        return (True);
    } else {
        return (False);
    }
}   /* feasible_move */

/* move the selected cell */
void move_cell(selected_cell_t scell[],
               allele tchrom[],
               cells_t cells[],
               parts_t parts[])
{
    parts[scell[0].from_part].pmax_cells--;
    parts[scell[0].to_part].pmax_cells++;
    parts[scell[0].from_part].pcurr_size -=
        cells[scell[0].mov_cell_no].cweight;
    parts[scell[0].to_part].pcurr_size +=
        cells[scell[0].mov_cell_no].cweight;
    tchrom[scell[0].mov_cell_no] = scell[0].to_part;
}   /* move_cell */

/* update gains after a move */
void update_gains(selected_cell_t scell[],
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  corn_t ncells[],
                  cells_info_t cells_info[],
                  allele tchrom[])
{
    int net_ptr = cells[scell[0].mov_cell_no].netlist;
    int mov_cell_no = scell[0].mov_cell_no;
    int from_part = scell[0].from_part;
    int to_part = scell[0].to_part;

    /* for each neighor net */
    for (int i = 0; i < cells[mov_cell_no].cno_nets; i++) {

        int net_no = cnets[net_ptr + i].corn_no;
        int net_weight = nets[net_no].nweight;
        int cell_ptr = nets[net_no].celllist;

        /* do operations before the move */
        if (nets[net_no].npartdeg[from_part] == nets[net_no].nno_cells) {

            update1(False, from_part, mov_cell_no, 
                    cell_ptr, net_no, net_weight,
                    nets, ncells, cells_info);

        } else if (nets[net_no].npartdeg[from_part] == (nets[net_no].nno_cells - 1)) {

            update2(False, from_part, mov_cell_no, 
                    cell_ptr, net_no, net_weight,
                    nets, ncells, cells_info, tchrom);
            
        }   /* else */

        /* update net info */
        nets[net_no].npartdeg[from_part]--;
        nets[net_no].npartdeg[to_part]++;

        /* do operations after the move */
        if (nets[net_no].npartdeg[to_part] == nets[net_no].nno_cells) {

            update1(True, to_part, mov_cell_no, 
                    cell_ptr, net_no, net_weight,
                    nets, ncells, cells_info);

        } else if (nets[net_no].npartdeg[to_part] == (nets[net_no].nno_cells - 1)) {

            update2(True, to_part, mov_cell_no, 
                    cell_ptr, net_no, net_weight,
                    nets, ncells, cells_info, tchrom);

        }   /* else */
    }   /* for i */
}   /* update_gains */

/* update gain of a cell */
void update1(int flag,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             nets_t nets[],
             corn_t ncells[],
             cells_info_t cells_info[])
{
    int i = 0;
    while (i < nets[net_no].nno_cells) {
        int other_cell = ncells[cell_ptr + i].corn_no;
        if (other_cell != mov_cell_no) {

            if (flag == False) {
                cells_info[other_cell].mgain[dest_part] -= net_weight;
            } else {
                cells_info[other_cell].mgain[dest_part] += net_weight;
            }

        }   /* if */

        i++;
    }   /* while */
}   /* update1 */

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
             allele tchrom[])
{
    int other_cell, other_part_no;

    int found = False;
    int i = 0;
    while ((i < nets[net_no].nno_cells) && (!found)) {

        other_cell = ncells[cell_ptr + i].corn_no;
        other_part_no = tchrom[other_cell];
        if ((other_cell != mov_cell_no) &&
            (other_part_no != dest_part)) {
            found = True;
        }   /* if */

        i++;
    }   /* while */

    if (!found) return;

    if (flag == False) {
        cells_info[other_cell].mgain[dest_part] -= net_weight;
    } else {
        cells_info[other_cell].mgain[dest_part] += net_weight;
    }
}   /* update2 */

/* copy pop structures */
void copy_pop(int nocells,
              int noparts,
              ind_t from_pop[], 
              ind_t to_pop[])
{
    for (int i = 0; i < nocells; i++) {
        to_pop[0].chrom[i] = from_pop[0].chrom[i]; 
    }   /* for i */
    for (int i = 0; i < noparts; i++) {
        to_pop[0].parts[i].pmax_cells = from_pop[0].parts[i].pmax_cells;
        to_pop[0].parts[i].pmax_size = from_pop[0].parts[i].pmax_size;
        to_pop[0].parts[i].pmin_size = from_pop[0].parts[i].pmin_size;
        to_pop[0].parts[i].pcurr_size = from_pop[0].parts[i].pcurr_size;
        to_pop[0].parts[i].pratio = from_pop[0].parts[i].pratio;
    }   /* for i */
    to_pop[0].incost = from_pop[0].incost;
}   /* copy_pop */

/* copy nets structures */
void copy_nets(int nonets,
               int noparts,
               nets_t from_nets[],
               nets_t to_nets[])
{
    for (int i = 0; i < nonets; i++) {
        to_nets[i].nno_cells = from_nets[i].nno_cells;
        to_nets[i].nweight = from_nets[i].nweight;
        to_nets[i].celllist = from_nets[i].celllist;
        for (int j = 0; j < noparts; j++) {
            to_nets[i].npartdeg[j] = from_nets[i].npartdeg[j];
        }        
    }
}  /* copy_nets */

/* EOF */
