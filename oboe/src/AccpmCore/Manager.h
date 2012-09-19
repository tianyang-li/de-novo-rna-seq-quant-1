// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef MANAGER_H
#define MANAGER_H

#include "AccpmDynMatrix.h"
#include "Parameters.h"
#include "ExitCode.h"

#include <map>
using std::map;

#ifdef SERIALIZATION
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#endif

/**
 *  @defgroup AccpmCore
 *  Group for managing the core Accpm functions 
 *  of creating the localization set and generating a query point.
 */
/**
 * @file Manager.h
 * @ingroup AccpmCore
 */

namespace Accpm 
{
 
/**
 * Manages the point generated - mirrors the PointGenS structure
 * Contains the decision variables as in PointGenS.DecVar 
 * @ingroup AccpmCore */  

  class PointGen {
  
  public :
    PointGen(const Parameters &param);
    ~PointGen() {};
    void setRho(double rho) { _rho = rho; }

    void updateATQA(SymmetricMatrix &H);

    AccpmVector _s;
    AccpmVector _sOld;
    AccpmVector _ds;
    AccpmVector _x;
    AccpmVector _xScaled;
    
    double _ws;// Ref GAC: omega
    double _zs; // Ref GAC: tau
    double _ss; // Ref: GAC: sigma
    AccpmVector _q;
    double _rho; // Ref: MATLAB code to replace PointGenS.rhoOld
    AccpmGenMatrix _ATQA; 

#ifdef SERIALIZATION

    friend class boost::serialization::access;
    template<class Archive> 
      void serialize(Archive &ar, const unsigned int file_version)
      {
	ar & _s;
	ar & _sOld;
	ar & _ds;
	ar & _x;
	ar & _xScaled;
	ar & _ws & _zs & _ss & _q & _rho;
      }
#endif
  };

   typedef std::pair<const AccpmVector *, int> AccpmVectorIntPair;
/**
 * Manages the data generated during Accpm iterations.
 *
 * @ingroup AccpmCore */  

  class Manager {
 
 
#ifdef SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> 
      void save(Archive &ar, const unsigned int file_version) const
      {
	std::cout << "Serializing Manager..." << std::endl;
	ar & _phase2 & _prevPhase2;
	ar & _rho & _objectiveFunction & _objLB & _objUB;
	ar & _relativeGap & _weightEpigraphCut;
	ar & _currentPointy & _bestPointy & _currentPointz;
	ar & _rhsCoef & _subProblemIndex & _cutOccurrence &_activeCutNorm;
	ar & _activeCuts;
	ar & _pointGen;
	ar & _sizeOfs & _numCuts & _hasSmoothOracle;
      }

    template<class Archive> 
      void load(Archive &ar, const unsigned int file_version)
      {
	std::cout << "Unserializing Manager..." << std::endl;
	ar & _phase2 & _prevPhase2;
	ar & _rho & _objectiveFunction & _objLB & _objUB;
	ar & _relativeGap & _weightEpigraphCut;
	ar & _currentPointy & _bestPointy & _currentPointz;
	ar & _rhsCoef & _subProblemIndex & _cutOccurrence & _activeCutNorm;
	ar & _activeCuts;
	ar & _pointGen;
	ar & _sizeOfs & _numCuts & _hasSmoothOracle;
	postprocess();
      }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    virtual void postprocess() {
      _activeCutsM = _activeCuts.getM();
      assert (_numCuts == _activeCutsM.size(1));
      assert (_numCuts == _subProblemIndex.size());
      assert (_numCuts == _cutOccurrence.size());
      assert (_numCuts == _rhsCoef.size());
      for (int i = 0; i < _numCuts; ++i) {
	AccpmVector *cut = new AccpmVector(_activeCutsM.getColumn(i));
	AccpmVectorIntPair *cutp = new AccpmVectorIntPair(cut,(int) _subProblemIndex(i));
	_cutSet[cutp] = ++_currentCutId;
      }
    }

#endif
     
  protected:
	
    struct ltAccpmVectorPair
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

      bool operator()(const AccpmVectorIntPair *v1, const AccpmVectorIntPair *v2) const
      {
		  //check the subproblem index
		  int v1s = v1->second;
		  int v2s = v2->second;
		  
		  if (v1s < v2s) {
	        return true;
		  }
		  if (v1s > v2s) {
		    return false;
		  }
		  v1s = v1->first->size();
		  v2s = v2->first->size();

		  if (v1s < v2s) {
	        return true;
		  }
		  if (v1s > v2s) {
		    return false;
		  }
		  return findFirstSmallerIndex(v1->first, v2->first);
      }
    };
   
    typedef map<const AccpmVectorIntPair *, int,  ltAccpmVectorPair> CutSet;  
    const Parameters *_param;

    bool _phase2;
    bool _prevPhase2;
    bool _currentPointIsFeasible;
    /** 
     * The main iteration counter for the query point generator.
     * At each of the outer iteration the Oracle::eval() is called.
     **/
    int _numOuterIteration; 
    
    /**
     * The total number of inner iterations done so far.
     * The inner iterations refer to the number of steps taken to solve the
     * Newton system of equations.
     * Typically this number is 2-3 per outer iteration.
     **/
    int _numInnerIteration; // It is the total number of inner iterations
    
    ExitCode _exitCode;
    double _rho;
    double _objectiveFunction;
    double _objLB;
    double _objUB;
    double _relativeGap;
    double _weightEpigraphCut;

    AccpmVector _currentPointy;
    AccpmVector _bestPointy;
    AccpmVector _currentPointz;
    
    AccpmVector _rhsCoef;
    AccpmVector _subProblemIndex;
    AccpmVector _cutOccurrence;
    AccpmVector  _activeCutNorm;
    AccpmDynMatrix _activeCuts;
    AccpmGenMatrix _activeCutsM;
    PointGen _pointGen;
    
    double _f2; //smoothObjFunctionVal;
    AccpmVector *_b; 
    AccpmVector *_df2; //smoothGrad; 
    AccpmGenMatrix *_d2f2; // smoothHessian;
    //bool _diagHessian;

    int _sizeOfs;
    int _numCuts;
    bool _hasSmoothOracle;

    /*
     * Hash for the cuts to help filtering duplicate cuts
     */
    CutSet _cutSet;
    int _currentCutId;

    virtual void init(void);
    virtual void update();
    virtual double convexityFix(const AccpmVector &functionValue,const AccpmGenMatrix &subGradients, 
			const AccpmVector &subProblemIndex,
			const AccpmVector &upperFunctionValues);
    virtual void updateRhs(const AccpmVector &functionValue, const AccpmGenMatrix &subGradients);
    virtual void updateRho(void);
    virtual void updateEpigraphWeight(void);
    
    virtual int addCut(const AccpmVectorIntPair *cut, int subProblemIndex, double weight, double rhs);
    int findCut(const AccpmVectorIntPair &v) const;
    virtual int processCuts(const AccpmVector &functionValue, const AccpmGenMatrix &newCuts, 
		    const AccpmVector &subProblemIndex, const AccpmVector &upperFunctionValues);
    virtual void updateQ();
    void updateRelativeGap();

  public:
    Manager(const Parameters *param);
    virtual ~Manager();
    virtual const AccpmVector& getCurrentY() const;
    virtual const AccpmVector& getActiveY() const;
    virtual int getActiveCuts(AccpmGenMatrix &cuts, bool activeVar = false) const;
    virtual int getActiveCuts(const AccpmGenMatrix *&cuts, bool activeVar = false) const;

    virtual const AccpmVector& getBestY(bool activeVar = false) const;
    inline const AccpmVector& getCurrentZ() const { return _currentPointz; }
    inline const AccpmVector& getRhsCoef() const { return _rhsCoef; }
    inline const AccpmVector& getSubProblemIndex() const { return _subProblemIndex; }
    inline const AccpmVector& getCutOccurrence() const { return _cutOccurrence; }
    inline const AccpmVector& getActiveCutNorm() const { return _activeCutNorm; }

    inline const AccpmVector& getCurrentX() const { return _pointGen._xScaled; }
           
    inline PointGen &getPointGen() { return _pointGen; }
    inline const PointGen &getPointGen() const { return _pointGen; }
    inline const AccpmGenMatrix &getATQA() const { return _pointGen._ATQA; }
    inline AccpmGenMatrix &getATQA() { return _pointGen._ATQA; }
    
    inline double getObjLB() const { return _objLB; }
    inline double getObjUB() const { return _objUB; }
    double getRelativeGap() const;

    inline double getWtEpigraphCut() const { return _weightEpigraphCut; }

    virtual double getSmoothObj(bool activeVar = false) const;
    virtual const AccpmVector *getSmoothGradient(bool activeVar = false) const;
    virtual const AccpmGenMatrix *getSmoothHessian(bool activeVar = false) const;
     
    inline bool inPhase2() const { return _phase2; }
    inline bool previousPhase2() const { return _prevPhase2; }
    inline void updatePreviousPhase2() { _prevPhase2 = _phase2; }
    inline bool isCurrentPointFeasible() const { return _currentPointIsFeasible; }
    inline ExitCode getExitCode() const { return _exitCode; }
    inline void setExitCode(ExitCode e) { _exitCode = e; }
   
    inline double getRho() const { return _rho; }
    inline int getSizeOfs() const { return _sizeOfs; }
    inline void setSizeOfs(int n) { _sizeOfs = n; }
    
    virtual int update1(const AccpmVector &functionValue, const AccpmGenMatrix &subGradients, 
		const AccpmGenMatrix &subProblemIndex,
		const AccpmVector *upperFunctionValues = 0);
    virtual int update2(double functionValue, const AccpmGenMatrix &subGradients, 
		const AccpmGenMatrix &hessian);
    
    double computeBallConstraint(const AccpmVector &y) const;
    void computeBallConstraint(const AccpmVector &y, AccpmVector &result) const;
    void computeBox1Constraint(const AccpmVector &y, AccpmVector &result) const;
    void computeBox2Constraint(const AccpmVector &y, AccpmVector &result) const;
    
    virtual void updateVariables(const AccpmVector &y, const AccpmVector &z, const AccpmVector &s,
				 double zs, double s0, double ss,
				 const AccpmVector &sOld, double sOld0,
				 const AccpmVector &ds, double ds0);
    virtual void updateY(const AccpmVector &y);

    int getNumOuterIteration() const { return _numOuterIteration; }

    void updateInnerIterations(int n) { _numInnerIteration += n; }
    int getNumInnerIterations() const { return _numInnerIteration; }

    virtual void updateLB(double lBound);

    virtual const AccpmVector *getVariableLB(bool activeVar = false) const;
    virtual const AccpmVector *getVariableUB(bool activeVar = false) const;
    virtual const AccpmVector *getB(bool activeVar = false) const; 
    virtual const AccpmVector *getCenterBall(bool activeVar = false) const; 
    virtual double getRadiusBall() const;

    virtual const AccpmVector *getQ(bool activeVar = false) const;
    virtual int computeSmoothComponent(bool hessian = true);
    int callSmoothOracle(const AccpmVector &y, double &f2, 
			 AccpmVector &df2, AccpmGenMatrix *d2f2) const;

    virtual void updateActiveVariables();
    virtual void updateActiveVariablesForF2();

    virtual int getNumCuts() const;
    
    virtual void output(std::ostream &os) const;
       
    void check();
  };
}

#endif
