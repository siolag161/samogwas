# samogwas project

# Overview
The `samogwas` project is organized structurally as follows: 

* dependencies: external resources needed to build the project:
   * BOOST project
   * Pro-BT library
* document: all related document 
* core: core source code of the project.
  * include: headers
     * clustering: related to clustering algorithms involving distance/similarity matrix  
       * `DBSCAN` 
       * `CAST`
       * `K-means`
     * distance: interfaces for distance/similarity matrices
     * em: interfaces and implementations of multiple Expectation-Maximization algorithms
     * fltm: implementation of the FLTM algorithm
     * statistics: statistics-related code:
       * statistical tests for association/correlation ( `chi-squared`, `fisher's`, `g-test`, etc..)
       * information theory related code
     * utils
       * helpful helpers code
   * src: source code 
* projects:
  *  stand-alone or projects that depend on `SAMOGWAS core` project
* sample: sample dataset for test purposes
* tests: unit tests for `core` and evolving projects
