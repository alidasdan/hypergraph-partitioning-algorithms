/***********************************************************************
** I give ALATURKA in an "as is" basis. I do not say or imply that    **
** it will be useful for whatever you want to do with it. You can     **
** use ALATURKA free of charge in academic research and teaching.     **
** For any commercial use, you must get a written permission from me. **
**                                                                    **
** (C) Ali Dasdan, 1991-                                              **
***********************************************************************/
/* FOR HYPERGRAPHS */
# define MAX_POP            1     /* max numbr of individuals in population */
# define NIL               -1     /* point to nowhere */
# define STR_SIZE          30     /* string size */

/* cells type */
typedef struct cells_st {
  int  cno_nets;           /* total number of nets on the cell */
  int  cno_inets;          /* number of internal nets */
  int  cweight;            /* weight of cell */
  int  netlist;            /* pointer to nets on the cell */
}   cells_type;
 
/* nets type */
typedef struct nets_st {
  int  nno_cells;          /* total number of cells on the net */
  int  nweight;            /* weight of node */
  int  celllist;           /* pointer to cells on the net */
  int  npartdeg [MAX_PARTS];  /* #cells of the net in each part */
}   nets_type;

/* nets info type */
typedef struct nets_info_st {
  int  npartdeg [MAX_PARTS];  /* #cells of the net in each part */
}   nets_info_type;
 
/* cell cornlist type - CorN = Cell or Net */
typedef struct corn_st {
  int corn_no;              /* cell or net number */
}   corn_type;

/* allele definition - for compatibility purposes */
typedef int allele;

/* partition type */
typedef struct parts_st {
  int  pmax_cells;         /* maximum number of cells in partition */
  int  pmax_size;        /* maximum size of part */
  int  pmin_size;        /* minimum size of part */
  int  pcurr_size;       /* current size of part */
  float pratio;            /* pmax_size / totsize */
}   parts_type;

/* individual type */
typedef struct pop_st {
  allele         chrom  [MAX_CELLS];   /* string holding partitions */
  parts_type     parts  [MAX_PARTS];   /* partition array */
  int            incost;               /* sum of net weights - cut cost */
}   ind_type;
 
/* temporary partition type */
typedef struct tparts_st {  /* used while creating partition celllists */
  char filled;              /* set if partition is filled with cells */
  int  pcells_inx;          /* an index to pcells */
}   tparts_type;

/* additional information for cells */
typedef struct cells_info_st {
  int  mgain [MAX_PARTS];   /* external costs of moving cell_no to all parts */ 
                            /* only mgain [part_no of cell_no] is internal cost */
}   cells_info_type; 

/* selected cell structure */
typedef struct selected_cell_st {
  int mov_cell_no;        /* current properties */
  int from_part;
  int to_part;
  int mov_gain;
}   selected_cell_type;

/* partition information type */
typedef struct parts_info_st {
  int  pmax_cells;         /* maximum number of cells in partition */
  int  pmax_size;        /* maximum size of part */
  int  pmin_size;        /* minimum size of part */
  int  pcurr_size;       /* current size of part */
}   parts_info_type;

/* precomputed exponential values */
typedef struct eval_st {
  double val;
}   eval_type;
