/**
 * @file   cond_ind_test.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Sun Mar 30 11:10:49 2008
 * 
 * @brief  R warp of cond_ind_test<>,
 *         Use "R CMD SHLIB cond_ind_test.cpp" to generate shared library
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <stdexcept>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "cond_ind_test.h"

extern "C" 
{

  /** 
   * wrap function of cond_ind_test<std::vector<int> >()
   * 
   * @param tb a real matrix with frequency as the last column
   * @param nvar number of different values as an integer vector 
   * @param i0 
   * @param i1 
   * 
   * @return p.value
   */
  SEXP cond_ind_test( SEXP tb, SEXP nvar, SEXP i0, SEXP i1 )
  {
    // get dimensions
    SEXP dim_attr = getAttrib(tb, R_DimSymbol );
    int m = INTEGER(dim_attr)[0];
    int n = INTEGER(dim_attr)[1]-1;
    double dev = 0.;
    double dev_dof = 0.;
    

    // wrap to stl form
    std::vector<std::pair<std::vector<int>, int> >
      stl_tb(m,std::pair<std::vector<int>, int>(std::vector<int>(n,0),1));

    for( int i=0; i<m; ++i )
      {
        for( int j=0; j<n; ++j )
          stl_tb[i].first[j] = INTEGER(tb)[i+m*j];
        stl_tb[i].second = INTEGER(tb)[i+m*n];
      }

    std::vector<int> stl_nvar(n);
    for( int j=0; j<n; ++j )
      stl_nvar[j] = INTEGER(nvar)[j];

    SEXP result = PROTECT(allocVector(REALSXP,3));
    REAL(result)[2] = cond_ind_test( stl_tb, stl_nvar, INTEGER(i0)[0] , INTEGER(i1)[0], dev, dev_dof );
    REAL(result)[0] = dev;
    REAL(result)[1] = dev_dof;
    UNPROTECT(1);

    return result;
  }
  
} // extern "C"
