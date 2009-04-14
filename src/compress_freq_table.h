/**
 * @file   compress_freq_table.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Fri Mar 28 19:14:02 2008
 * 
 * @brief C++ function to compress frequency table
 * 
 * 
 */

#ifndef COMPRESS_FREQ_TABLE_HPP
#define COMPRESS_FREQ_TABLE_HPP

#include <vector>
#include <algorithm>

/** 
 * Compress the frequency table
 * 
 * @param v original frequency table
 */
template<typename T>
void compress_freq_table( std::vector<std::pair<T, int> >& v )
{
  std::sort( v.begin(), v.end() );

  typedef typename std::vector<std::pair<T,int> >::iterator iter_t;
  
  iter_t last_v = v.begin()-1;  // pointer to the last unique entry
  iter_t cur_v  = v.begin();    // pointer to the current entry

  // compress the frequency table
  while( cur_v != v.end() )
    {
      if( cur_v->second != 0 )  
        {
          // create a new entry
          if( last_v == v.begin()-1 || cur_v->first != last_v->first )
            {
              last_v++;
              *last_v = *cur_v;
            }
          // update an existing entry
          else
            {
              last_v->second += cur_v->second;
            }
        }

      cur_v++;
    }

  v.erase( last_v+1, v.end() );
}

#endif  // COMPRESS_FREQ_TABLE_HPP

