# core source code

## Overview directory map

+  clustering: interfaces and implementation for clustering algorithms and partition classes
    +  clustering.hpp: clustering algorithm interfaces
    +  partition.hpp: definition of partition/clustering and clusters
    +  dbscan.hpp: interface for DBSCAN algorithm 
    +  cast.hpp: interface for CAST algorithm

+  distance:
    +  comparable.hpp: common interface for similarity/dissimilarity matrix interfaces
    +  dissimilarity.hpp: common interface for dissimilarity matrix functor
    +  information_dissimilarity.hpp: dissimilarity based on information theory (correlation, mutual information etc..)
    +  information_similarity.hpp: similarity based on information theory
    +  similarity.hpp: common interface for dissimilarity matrix functor 

+  em:
    +  core_em.hpp: common methods for EM algorithms   
    +  em_helper.hpp: helper methods for EM algorithms
    +  em.hpp: common interfaces for EM algorithms  
    +  naive_bayes_em.hpp: implementation of the EM algorithm for naive-bayes model
+  fltm
+  statistics: 
+  utils:
    +  csv_parser.hpp: a simple tool for parsing csv-like data files. It is implemented in order to be used in a iterator-like manner.  
    +  matrix_utils.hpp: several methods for handing matrix objects
    +  type_utils.hpp: helpers for dealing with type and structures
