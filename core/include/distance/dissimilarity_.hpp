/****************************************************************************************
 * File: Diss.hpp // CS This is not the name of the file.
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com) // CS Affiliation
 * @date: 09/07/2014

 ***************************************************************************************/
#ifndef FLTM_DISS_HPP // CS Heterogeneous name assignement rule
#define FLTM_DISS_HPP

#include <vector>

#include "mut_info_distance.hpp" // CS How can this code work. This file does not even exist.
namespace fltm // CS Heterogeneous name assignement rule
{

struct Comp {
  virtual double operator()( const size_t varA, const size_t varB ) { return compute(varA, varB); }
  virtual double compute( const size_t varA, const size_t varB ) = 0;
  virtual size_t size() const = 0;
  virtual void invalidCache() = 0;
  // virtual void reset( std::vector< std::vector<int> >& dm, std::vector<int>& pos ) = 0;
};

struct Distance: public Comp {};


struct MutInfoDiss: public Distance {
  MutInfoDiss( std::vector< std::vector<int> >& dm, std::vector<int>& pos, unsigned maxPos, double thres ):
      m_diss(dm, pos, maxPos, thres) {}
  virtual double compute( const size_t varA, const size_t varB ) {
    return m_diss(varA,varB);
  }

  virtual size_t size() const { return m_diss.size(); }
  virtual void invalidCache() { m_diss.invalidCache(); }
 
private:
  MutInfoDistance< std::vector< std::vector<int> > > m_diss; // CS I expect the type MutInfoDistance to be defined in
                                                             // the unique included file "mut_info_distance.hpp"
                                                             // CS !!!!!!!!!!!!!!!!!!!!
                                                             // But I do not even see that the file mut_info_distance.hpp exists.
};


} // namespace fltmends here. fltm

/****************************************************************************************/
#endif // FLTM_DISS_HPP
