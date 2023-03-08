/*
  ---------------------------------------------------------------------
  This file is part of GRIP (Row Inner-Product Graph partitioner for ABCD)
  Copyright (c) 2018-   Fahreddin Sukru TORUN
  ---------------------------------------------------------------------
  For license info, please see the README file in the main directory.
  ---------------------------------------------------------------------

*/
#include <abcd.h>
#include "mat_utils.h"
#include <iostream>
#include <fstream>
#include <vect_utils.h>
#include <queue>
#include <vector>
#include <metis.h>

// Column Sparsing factor for dropping elements
// 0: no Sparsifying
// x: sqrt(n) * x largest elements are kept
double ColSparsing_thresh=0;

// scale factor for casting integer
int metis_scaleMax = 100;

//seed for partitioning
int metis_seed = -1;

// for dropping some tiny values in graph after AAT
double metis_epsilon = 0; // 1e-16;

// for calling METIS_PartGraphRecursive (optional)
int metis_recursive = 0; 

CompRow_Mat_double ColumnSparsifying(CompRow_Mat_double &inMtx);

void GRIP_Partitioning(CompRow_Mat_double &Ddrop,idx_t nParts, idx_t* part, double imb);

