// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef ACCPM_VECTOR_H
#define ACCPM_VECTOR_H
#include "AccpmDefs.h"
#include <vector>
#include "lavd.h"

typedef double Real;
typedef LaVectorDouble RealVector;
typedef std::vector<Real> StdRealVector;

#ifdef SERIALIZATION
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#endif

/**
 *  @defgroup AccpmLA
 *  Group for managing the Linear Algebra classes and operations.
 */
/**
 * @file AccpmVector.h
 * @ingroup AccpmLA
 */
namespace Accpm 
{

/**
 * Class for handling vectors which are compatible
 * with BLAS and LAPACK.
 * @ingroup AccpmLA */  

  class AccpmVector : public RealVector {

#ifdef SERIALIZATION  
    friend class boost::serialization::access;
    template<class Archive> 
      void save(Archive &ar, const unsigned int file_version) const
      {
	int n = size();
	ar & n;
	for (int i = 0; i < size(); ++i) {
	  ar & *(addr() + i);
	}
      } 

    template<class Archive> 
      void load(Archive &ar, const unsigned int file_version)
      {
	int n;
	ar & n;
	resize(n, 1);
	for (int i = 0; i < n; ++i) {
	  ar & *(addr() + i);
	}
      } 
    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  public:
    AccpmVector();
    AccpmVector(int n);
    AccpmVector(double* v, int n);
    AccpmVector(const LaGenMatDouble &s);
    virtual ~AccpmVector() {};

    AccpmVector &operator=(double s);
    AccpmVector& operator=(const LaGenMatDouble &s);
    bool operator==(const AccpmVector &v) const;
    void copy(const StdRealVector &v1); 
    void append(const AccpmVector &v);
    void append(const double entry); 
    int deleteElem(int id);
    void negate();
    void invert();
    double min(int *index = 0) const;
    double max(int *index = 0) const;
    /** 
     * Element by element mutliplication .* or times operation of MATLAB
     */
    void times(const AccpmVector &v); 
    /** 
     * Element by element division ./ or rdivide operation of MATLAB
     */
    void rdivide(const AccpmVector &v);
    double sum() const;
  };

 struct ltvector
    {
      bool findFirstSmallerIndex(const AccpmVector *v1, const AccpmVector *v2) const
      {
	for (int i = 0; i < v1->size(); ++i) {
	  if (DBL_LT((*v1)(i), (*v2)(i))) {
	    return true;
	  }
	  if (DBL_LT((*v2)(i), (*v1)(i))) {
	    return false;
	  } 
	}
	return false;
      }

      bool operator()(const AccpmVector *v1, const AccpmVector *v2) const
      {
	int v1s = v1->size();
	int v2s = v2->size();

	if (v1s < v2s) {
	  return true;
	}
	if (v1s > v2s) {
	  return false;
	}
/*
	double v1Value = v1->sum();
	double v2Value = v2->sum();
	if (DBL_LT(v1Value, v2Value)) {
	  return true;
	}
	if (DBL_GT(v1Value, v2Value)) {
	  return false;
	}
*/
	return findFirstSmallerIndex(v1, v2);
      }
    };
}

#endif
