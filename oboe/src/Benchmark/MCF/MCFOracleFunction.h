// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef MCF_ORACLE_FUNCTION_H
#define MCF_ORACLE_FUNCTION_H

#include "Oracle.h"
#include "AccpmGenMatrix.h"
#include <string>
#include <map>
#include <vector>

using std::string;
using std::map;
using std::vector;
using namespace Accpm;

class Arc {

public:
  unsigned int _from; // node id of from node
  unsigned int _to;    // node id of to node
  double _capacity;
  double _cost;
  Arc(unsigned int from, unsigned int to, double capacity, double cost) : 
    _from(from), _to(to), _capacity(capacity), _cost(cost) {}
};

class MCFData {

 public:
  typedef map<int, Arc *> ArcIdMap;
  typedef vector<int> IntVector;
  enum FunctionType { LINEAR = 'l', KLEINROCK = 'k', UNDEFINED = 'U'};
  
  unsigned int _numNodes;
  unsigned int _numArcs;
  unsigned int _numCommodities;
  unsigned int _numSourceNodes;
   
  AccpmGenMatrix _C;
  AccpmGenMatrix _D;
  AccpmVector _capacity;
  AccpmVector _cost;
  AccpmVector _demand; // Demands from third column of D matrix == D(:,3)

  FunctionType _type;
  int _stop; // signal early termination

  ArcIdMap _carcIdMap; // arcs are ordered in order of from node and then to node
  IntVector _clistAdj; 
  // List of to nodes in order of arcs; 
  IntVector _cindexAdj; 
  // For each from node, in order of its id,  gives the number of first arc not incident on it 
  
  ArcIdMap _darcIdMap; // arcs are ordered in order of from node and then to node for the demands 
  IntVector _dlistAdj; // Index Adjacency for demands == D(:,2)
  IntVector _dindexAdj; //Adjacency list for demands == IndCom from Matlab code
  
  AccpmVector _sol; // Keeps current value of Ax
  double _sumX; // Sum of x needed to scale if sum != 1
 private:
  int readMatrix(const char *fileName, AccpmGenMatrix &matrix);
  void buildInternalData();
  void buildInternalNetworkData();
  void buildInternalDemandData(); 
 
 public:
  MCFData(const char *networkFile, const char *dataFile, int numArcs, int numCommodities, char type);
  unsigned int getNumSources() const { return _numSourceNodes; }

  virtual ~MCFData();
};

class MCFOracleFunction : public OracleFunction {
 private:
  int _aggreg;
  MCFData *_data;

  int dijkstra(const double *value, double *grad);
  
public:

  MCFOracleFunction();
  MCFOracleFunction(const char *networkFile, const char *demandFile, int numArcs, int numCommodities,
		    char type);
  virtual ~MCFOracleFunction();
    
  int getNumNodes() const { return _data->_numNodes; }
  const AccpmVector &getCapacity() const;
  const MCFData *getMCFData();

  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info);
  int updateAx(const AccpmVector &sol, double sum = 1);
  void setAggregationLevel(int level);

};

#endif
