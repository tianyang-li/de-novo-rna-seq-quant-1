// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ACCPM_DYNMATRIX_H
#define ACCPM_DYNMATRIX_H

#include "AccpmGenMatrix.h"

#ifdef SERIALIZATION
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#endif

/**
 * @file AccpmGenMatrix.h
 * @ingroup AccpmLA
 */

namespace Accpm 
{

/**
 * Class for handling dynamic matrices which are capable of resizing
 * @ingroup AccpmLA */  

  class AccpmDynMatrix {
  private:
    bool _preallocateRows;
    bool _preallocateCols;
    int _currentRow;
    int _currentCol;

    AccpmGenMatrix _matrix;

    /** resize the matrix keeping as much of the matrix as possible
     */
    AccpmDynMatrix& resize(int m, int n);
#ifdef SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> 
      void serialize(Archive &ar, const unsigned int file_version)
      {
	ar & _preallocateRows & _preallocateCols & _currentRow & _currentCol;
	ar & _matrix;
      } 
#endif
  public:
    // Operations of RealMatrix
    AccpmDynMatrix(bool preallocateRows = true, bool preallocateCols = true, int m = 0, int n = 0); 
    AccpmDynMatrix(const AccpmGenMatrix &rhs, bool preallocateRows = true, bool preallocateCols = true);
    AccpmDynMatrix(const AccpmDynMatrix &);
   
    AccpmDynMatrix& operator=(double s);
    AccpmDynMatrix& operator=(const AccpmDynMatrix &);

    virtual ~AccpmDynMatrix();
    
    /**
     * return the currently used matrix
     */
    AccpmGenMatrix getM() const; 

    /** adds a column in the end and returns its id in the Matrix
     */
    int addColumn(const AccpmVector &v); 
    /** adds a row in the end and returns its id in the Matrix
     */
    int addRow(const AccpmVector &v);
    int deleteColumn(int id);

    int size(int n) const;
    std::ostream& Info(std::ostream& os) { _matrix.Info(os); return os; }
    friend std::ostream& operator<<(std::ostream& os, const AccpmDynMatrix &mat) { os << mat._matrix; return os; }
  };
}

#endif
