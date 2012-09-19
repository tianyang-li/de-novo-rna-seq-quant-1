// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ACCPM_GENMATRIX_H
#define ACCPM_GENMATRIX_H

#include "gmd.h"
#include "symd.h"

#include "AccpmVector.h"
#include "config.h"
#ifdef SERIALIZATION
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#endif

/**
 * @file AccpmGenMatrix.h
 * @ingroup AccpmLA
 */

//typedef LaSymmMatDouble SymmetricMatrix;

namespace Accpm 
{
  typedef LaGenMatDouble RealMatrix;
  typedef LaSymmMatDouble SymmetricMatrix;
/**
 * Class for handling matrices which are compatible
 * with LAPACK++ */  

  class AccpmGenMatrix : public RealMatrix {

#ifdef SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> 
      void save(Archive &ar, const unsigned int file_version) const
      {
	ar & size(0) & size(1);
	for (int i = 0; i < size(0) * size(1); ++i) {
	  ar & *(addr() + i);
	}
      }
    
    template<class Archive> 
      void load(Archive &ar, const unsigned int file_version) 
      {
	int m, n;
	ar & m & n;
	resize(m, n);
	for (int i = 0; i < m * n; ++i) {
	  ar & *(addr() + i);
	}
      }
    BOOST_SERIALIZATION_SPLIT_MEMBER() 
#endif
  
   public:
    // Operations of RealMatrix
    AccpmGenMatrix();
    AccpmGenMatrix(int m, int n);
    AccpmGenMatrix(double *v, int m, int n, 
		   bool row_ordering);
    AccpmGenMatrix(const AccpmGenMatrix &);
    AccpmGenMatrix(const RealMatrix &);
    AccpmGenMatrix(const AccpmVector &);

    virtual ~AccpmGenMatrix() {};

#if (defined(LINUX)) // The following do not work correctly on Windows .NET 2003
    inline AccpmGenMatrix operator()(const LaIndex& I, const LaIndex& J) const  
	{ return RealMatrix::operator()(I,J); }
    inline AccpmGenMatrix operator()(const LaIndex& I, const LaIndex& J) 
	{ return RealMatrix::operator()(I,J); }
    inline double& operator()(int i, int j) { return RealMatrix::operator()(i,j); }
    inline double& operator()(int i, int j) const { return RealMatrix::operator()(i,j); }                             
#endif

	//New AccpmGenMatrix operations
    AccpmGenMatrix* transpose() const;

    AccpmGenMatrix& operator=(double s);
    AccpmGenMatrix& operator=(const AccpmGenMatrix &s);
    RealVector getColumn(int i) const;
    RealVector getRow(int i) const;

    void scaleColumn(int i, double d);
    /**
     * Scale column i by a vector 
     * col(i) .* d
     */
    void scaleColumn(int i, const AccpmVector &d);
    void assignColumn(int colId, const AccpmVector &v);
    void assignRow(int rowId, const AccpmVector &v);
    /**
     * Equivalent of Blas_Add_Mult for Matrices 
     */
    void addMult(double scale, const AccpmGenMatrix &b);
    /**
     * Equivalent of Blas_Scale for Matrices 
     */
    void scale(double scale);
  };
}

#endif
