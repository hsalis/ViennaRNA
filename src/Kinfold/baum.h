/*
  Last changed Time-stamp: <2001-08-02 14:50:41 xtof>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: baum.h,v 1.1 2001/08/02 16:48:58 xtof Exp $
*/

#ifndef BAUM_H
#define BAUM_H

/* used in main.c */
extern void ini_or_reset_rl(void);
extern void move_it(void);
extern void clean_up_rl(void);
extern void kinfold_rebuild_current_state(void);
extern int  kinfold_can_pair_positions(int i, int j);
extern int  kinfold_can_pair_positions_raw(int i, int j);
extern int  kinfold_is_unpaired_position(int i);
extern int  kinfold_eval_structure_dcal(const char *structure);

/* used in nachbar.c */
extern void update_tree(int i,int j);

#endif
