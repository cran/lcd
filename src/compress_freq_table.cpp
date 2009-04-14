/**
 * @file   compress_freq_table.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Sun Mar 30 10:00:02 2008
 * 
 * @brief  R Use "R CMD SHLIB compress_freq_table.cpp" to generate shared library
 * 
 * 
 */

#include <vector>
#include <stdexcept>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "compress_freq_table.h"

extern "C" 
{

  /** 
   * wrap function of compress_freq_table<int>() 
   * 
   * @param tb an integer matrix with frequency as the last column
   * 
   * @return compressed freq table as an integer matrix
   */
  SEXP compress_freq_table( SEXP tb )
  {
    // get dimensions
    SEXP dim_attr = getAttrib(tb, R_DimSymbol );
    int m = INTEGER(dim_attr)[0];
    int n = INTEGER(dim_attr)[1]-1;

    // wrap to stl form
    std::vector<std::pair<std::vector<int>, int> >
      stl_tb(m,std::pair<std::vector<int>, int>(std::vector<int>(n,0),1));

    for( int i=0; i<m; ++i )
      {
        for( int j=0; j<n; ++j )
          stl_tb[i].first[j] = INTEGER(tb)[i+m*j];
        stl_tb[i].second = INTEGER(tb)[i+m*n];
      }

    // compress the frequency table
    compress_freq_table( stl_tb );

    // wrap back to SEXP
    m = stl_tb.size();
    
    SEXP c_tb = PROTECT(allocMatrix(INTSXP, m, n+1 ));
    for( int i=0; i<m; ++i )
      {
        for( int j=0; j<n; ++j )
          INTEGER(c_tb)[i+m*j] = stl_tb[i].first[j];
        INTEGER(c_tb)[i+m*n] = stl_tb[i].second;
      }
    
    UNPROTECT(1);

    return c_tb;
  }

} // extern "C"
