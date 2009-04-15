/**
 * @file   cond_ind_test.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Fri Mar 28 21:04:52 2008
 * 
 * @brief  Conditional independence test
 * 
 * 
 */

#ifndef COND_IND_TEST_HPP
#define COND_IND_TEST_HPP

#include <vector>
#include <map>

#include <cmath>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/** 
 * Conditional independence test
 *
 * M is obtained from compress_freq_table()
 * 
 * @param M measurement matrix with frequency, data must be the same type.
 * @param nvar number of different variables of each column
 * @param i index of first variable 
 * @param j index of second variable
 *
 * @return p.value
 * 
 */
template<typename T>
double cond_ind_test( const std::vector<std::pair<std::vector<T>, int> >& M,
                      const std::vector<int>& nvar,
                      int i, int j, double & dev, double & dev_dof)
{
  const int m = M.size();
  const int n = M[0].first.size();

  // build indices
  
  std::map<T, int> ii_i;  // i-th -> int
  std::map<T, int> ii_j;  // j-th -> int
  std::map<std::vector<T>, int> ii_k; // K-th -> int
  std::vector<std::vector<int> > ii(m);     // index of each measurement

  std::vector<int> max_ijk(3,0); // number of different values

  std::vector<T> x_k(n-2);
  
  for( int l=0; l<m; ++l )
  {
    const std::vector<T>& x = M[l].first;

    int t = 0;
    for( int s=0; s<n; ++s )
      if( s != i && s != j )
        x_k[t++] = x[s];
    
    if( !ii_i.count( x[i] ) )
      ii_i[x[i]] = (max_ijk[0]++);
    if( !ii_j.count( x[j] ) )
      ii_j[x[j]] = (max_ijk[1]++);
    if( !ii_k.count( x_k  ) )
      ii_k[x_k ] = (max_ijk[2]++);
    
    std::vector<int> cur_ijk(3);
    cur_ijk[0] = ii_i[x[i]];
    cur_ijk[1] = ii_j[x[j]];
    cur_ijk[2] = ii_k[x_k ];

    ii[l] = cur_ijk;
  }

  // count appearance

  std::map<std::vector<int>, int> f_ik; // N(i,*,k)
  std::map<std::vector<int>, int> f_jk; // N(*,j,k)
  std::vector<int> f_k(max_ijk[2]);           // N(*,*,k)

  for( int l=0; l<m; ++l )
    {
      int k = ii[l][2];
      std::vector<int> ik(2), jk(2);
      ik[0] = ii[l][0];
      jk[0] = ii[l][1];
      ik[1] = jk[1] = k;
  
      f_k[k]   += M[l].second;
      f_ik[ik] += M[l].second;
      f_jk[jk] += M[l].second;
    }

  // conditional independence test
  
  dev = 0.;
  for( int l=0; l<m; ++l )
    {
      int k = ii[l][2];
      std::vector<int> ik(2), jk(2);
      ik[0] = ii[l][0];
      jk[0] = ii[l][1];
      ik[1] = jk[1] = k;
      
      dev += M[l].second * log( 1. * M[l].second*f_k[k] / ( f_ik[ik]*f_jk[jk] ) );
    }

  dev *= 2.;

  double prod = 1.;
  for( int l=0; l<n; ++l )
    prod *= nvar[l];

  dev_dof = prod - prod/nvar[i] - prod/nvar[j] + prod/nvar[i]/nvar[j];
    
  // call R to compute p_value = 1-pchisq(dev,devdf)  
  return 1.-pchisq(dev, dev_dof, 1, 0 );
}

#endif // COND_IND_TEST_HPP
