/*
  Last changed Time-stamp: <2010-06-24 14:50:01 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: main.c,v 1.5 2008/08/28 09:40:55 ivo Exp $
*/

#include "config.h"
#include <stdio.h>

#include "kinfold_api.h"

int
main(int argc, char *argv[])
{
  return kinfold_run_with_args(argc, argv, stdin, NULL, 0);
}
