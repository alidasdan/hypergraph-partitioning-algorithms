#ifndef AD_READINPUT_INCLUDED
#define AD_READINPUT_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* initialize cells array */
void init_cells(int nocells, cells_t cells[]);

/* initialize nets array */
void init_nets(int nonets,
               int noparts,
               nets_t nets[]);

/* initialize netlist pointers */
void init_netlist(int nonets,
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  corn_t ncells[]);

void read_hgraph_size(char fname[],
                      int  *nocells,
                      int  *nets,
                      int  *nopins);

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
                 corn_t  ncells[]);

#endif
