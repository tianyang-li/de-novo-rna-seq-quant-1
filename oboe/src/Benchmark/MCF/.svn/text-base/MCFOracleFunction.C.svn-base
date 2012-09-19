// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "MCFOracleFunction.h"
#include "AccpmBlasInterface.h"
#include <iostream>
#include <fstream>
#include <sstream> 
 
MCFData::MCFData(const char *networkFile, const char *demandFile, 
		 int numArcs, int numCommodities, char type) :
  _numNodes(0), _numArcs(numArcs), 
  _numCommodities(numCommodities), _numSourceNodes(0),
  _C(_numArcs, 4), _D(_numCommodities, 3), 
  _capacity(_numArcs), _cost(_numArcs), _demand(_numCommodities), 
  _type(LINEAR), _stop(0), _sumX(1)
{
  if (tolower(type) == 'k') {
    _type = KLEINROCK;
  }
  int status = 0;
  status |= readMatrix(networkFile, _C);
  status |= readMatrix(demandFile, _D);
  if (status) {
    exit(1);
  }
  buildInternalData();
}
MCFData::~MCFData()
{
  ArcIdMap::iterator iter;
  for (iter = _carcIdMap.begin(); iter != _carcIdMap.end(); ++iter) {
    Arc *arc = (*iter).second;
    delete arc;
  }
  for (iter = _darcIdMap.begin(); iter != _darcIdMap.end(); ++iter) {
    Arc *arc = (*iter).second;
    delete arc;
  }
}

int 
MCFData::readMatrix(const char *fileName, AccpmGenMatrix &matrix) 
{
  //std::cout << "Reading File:" << fileName << std::endl;
  int status = 0;
  std::ifstream fin(fileName);
  string line;
  std::istringstream instream;
  double val;
  if (fin.is_open()) {
    for (int i = 0; i < matrix.size(0); ++i) {
      if (getline(fin, line)) {
	instream.clear();
	instream.str(line); 
	for (int j = 0; j < matrix.size(1); ++j) { 
	  if ((instream >> val >> std::ws)) {
	    matrix(i, j) = val; 
	  } else {
	    std::cerr << "Error reading matrix from file: " << fileName 
		      << " at line:" << i+1 << " column: " << j + 1 << std::endl;
	    status = 1;
	  }
	}
	getline(fin, line); // eat up empty line which is there in current format
      } else {
	std::cerr << "Error reading matrix from file: " << fileName << std::endl;
	status = 1;
      }
    }
    } else {
    AccpmError("Cannot open file: ");
    std::cout << fileName << std::endl;
    status = 1;
  }
  return status;
}

void
MCFData:: buildInternalData() 
{
  buildInternalNetworkData();
  buildInternalDemandData();
  assert(_numArcs == _carcIdMap.size());
  assert(_numCommodities == _darcIdMap.size());
  assert(_numNodes + 1 == _cindexAdj.size());
  assert(_numNodes + 1 == _dindexAdj.size());

}
  
void
MCFData::buildInternalNetworkData() 
{
  AccpmVector c0 = _C.getColumn(0);
  AccpmVector c1 = _C.getColumn(1);
  AccpmVector c2 = _C.getColumn(2);
  AccpmVector c3 = _C.getColumn(3);
  
  _numNodes = (int)std::max(c0.max(), c1.max());
  // build the arc id Map
  for (unsigned int i = 0; i < _numArcs; ++i) {
    int from = int(c0(i) - 1);
    int to = (int)(c1(i) - 1);
    int arcId = from * _numNodes + to;
    if (_carcIdMap.find(arcId) == _carcIdMap.end()) {
      double capacity = c2(i);
      double cost = c3(i);
      _carcIdMap[arcId] = new Arc(from, to, capacity, cost);
    } else {
      std::cerr << "Network file error: Duplicate arc found: (" << from + 1 
		<< " , " << to + 1 << " ). Keeping the orignal." << std::endl;
    } 
  }
  ArcIdMap::iterator iter;
  unsigned int j = 0;
  for (iter = _carcIdMap.begin(); iter != _carcIdMap.end(); ++iter, ++j) {
    Arc *arc = (*iter).second;
    _clistAdj.push_back(arc->_to);
    _capacity(j) = arc->_capacity;
    _cost(j) = arc->_cost;
  }

  iter = _carcIdMap.begin();
  j = 0;
  for (unsigned int i = 0; i < _numNodes; ++i) {
    while (iter != _carcIdMap.end() && 
	   (*iter).second->_from + 1 == i) {
      ++iter;
      ++j;
    }
    _cindexAdj.push_back(j);
  }
  _cindexAdj.push_back(_numArcs);
}

void
MCFData::buildInternalDemandData() 
{
  AccpmVector d0 = _D.getColumn(0);
  AccpmVector d1 = _D.getColumn(1);
  AccpmVector d2 = _D.getColumn(2);
  
  // build the arc id Map
  for (unsigned int i = 0; i < _numCommodities; ++i) {
    int from = int(d0(i) - 1);
    int to = (int)(d1(i) - 1);
    int arcId = from * _numNodes + to;
      if (_darcIdMap.find(arcId) == _darcIdMap.end()) {
	double demand = d2(i);
	_darcIdMap[arcId] = new Arc(from, to, demand, 0);
      } else {
	std::cerr << "Demand file error: Duplicate arc found: (" << from + 1 
		  << " , " << to + 1 << " ). Keeping the orignal." << std::endl;
      } 
  }
  
  unsigned int j = 0;
  ArcIdMap::iterator iter;
  for (iter = _darcIdMap.begin(); iter != _darcIdMap.end(); ++iter, ++j) {
    Arc *arc = (*iter).second;
    _dlistAdj.push_back(arc->_to);
    _demand(j) = arc->_capacity;
  }
  
  iter = _darcIdMap.begin();
  j = 0;
  for (unsigned int i = 0; i < _numNodes; ++i) {
    while (iter != _darcIdMap.end() && 
	   (*iter).second->_from + 1 == i) {
      ++iter;
      ++j;
    }
    _dindexAdj.push_back(j);
  }
  _dindexAdj.push_back(_numCommodities);
  
  iter = _darcIdMap.begin();
  while(iter != _darcIdMap.end()) {
    unsigned int from = (*iter).second->_from;
    ++_numSourceNodes;
    while (iter != _darcIdMap.end() && from == (*iter).second->_from) {
      ++iter;
    }
  }
  // std::cout << "nbSourceNodes = " << _numSourceNodes << std::endl;
}

MCFOracleFunction::MCFOracleFunction() : OracleFunction() 
{

}

MCFOracleFunction::MCFOracleFunction(const char *networkFile, const char *demandFile, int numArcs, 
				     int numCommodities, char type) : 
  OracleFunction(), _aggreg(1)
{
  _data = new MCFData(networkFile, demandFile, numArcs, numCommodities, type);
}

MCFOracleFunction::~MCFOracleFunction()
{
  if (_data) {
    delete _data;
  }
}

const MCFData *
MCFOracleFunction::getMCFData()
{
  return _data;
}

/**
 * Based on Matlab mex implementation in versionM/Benchmark/MCF/Dijkstra.c
 * The difference is that the gradients are already negated and dont have to 
 * handled by the calling function.
 */
int
MCFOracleFunction::dijkstra(const double *value, double *grad)
{
  unsigned int n = _data->_numNodes;
  double *demandV = _data->_demand.addr();

  // Temporaries for Dijkstra
  unsigned int *Pred = new unsigned int[n];
  short int *StateNode = new short int[n];
  double *D = new double[n];
  int *IndexArcPred = new int[n];
  short int *ListNode = new short int[n];
  unsigned int *IndexListNode = new unsigned int[n];
  int CurrentSourceNode = 0;
  for(unsigned int i = 0; i < n; ++i) { 
    IndexArcPred[i] = 0;
  }
  /* On traite chacun des noeuds sources s */
  for(unsigned int s = 0; s < n; s++) { 
    /* cherche l'indice min et l'inidce max ou il y a une demande a partir du
       noeud source s */
    unsigned int minIndCom = _data->_dindexAdj[s];
    unsigned int maxIndCom = _data->_dindexAdj[s + 1];
    unsigned int nIndCom = maxIndCom - minIndCom; /* Nombre de demande a partir de s */
    
    /* si s est un noeud source, on cherche les ppc */
    if (nIndCom > 0) {
      /*Initialization pre process*/
      /* Initialization of the vector of distance D, state node and predecesseur  */
      for(unsigned int i = 0; i < n; i++) {
	D[i] = 1e20;
	StateNode[i] = 0;
	//Pred[i] = s + 1;
	Pred[i] = s;
      }
      /* */
      int MaxNode = 0;
      int min = _data->_cindexAdj[s];
      int max = _data->_cindexAdj[s + 1];
      for(int i = min; i < max; i++) {
	int k = _data->_clistAdj[i];
	D[k] = value[i];
	IndexArcPred[k] = i;
	StateNode[k] = 1;
	ListNode[MaxNode] = k;
	IndexListNode[MaxNode] = MaxNode + 1;
	MaxNode = MaxNode + 1;
      }
      IndexListNode[MaxNode - 1] = n;
      
      for(unsigned int j = 0; j < n; j++) {
	/* Find minimum distances Dmin in vector D */
	double Dmin = 1e15;
	int IndexPred = 0;
	unsigned int Index = 0;
	/*On verifie que le premier element de la liste n'a pas deja ete traite */
	int k = ListNode[0];
	if (StateNode[k] == 2)
	  Index = IndexListNode[Index];
	
	unsigned int IndexMinPred = 0;
	int IndexMin = 0;
	/* Recherche du noeud minimum */
	while(Index != n) {
	  k = ListNode[Index];
	  if(D[k] < Dmin) {
	    Dmin = D[k];
	    IndexMinPred = IndexPred;
	    IndexMin = Index;
	  }
	  IndexPred = Index;
	  Index = IndexListNode[Index];
	}
	
	if(Dmin != 1e15) {
	  /* Update D, ListNode and IndexListNode */
	  
	  /* Retire le noeud traite */
	  int k = ListNode[IndexMin];
	  StateNode[k] = 2;
	  /* Ajouter les noeuds adjacents au noeud traite */
	  int min = _data->_cindexAdj[k];
	  int max = _data->_cindexAdj[k + 1];
	  for(int i = min; i < max; i++) {
	    unsigned int adjacent = _data->_clistAdj[i];
	    if(adjacent != s) {
	      if (StateNode[adjacent] == 0) {
		StateNode[adjacent] = 1;
		ListNode[MaxNode] = adjacent;
		IndexListNode[MaxNode] = n;
		IndexListNode[MaxNode - 1] = MaxNode;
		MaxNode = MaxNode + 1;
		D[adjacent] = value[i] + D[k];
		IndexArcPred[adjacent] = i;
		//Pred[adjacent] = k + 1;
		Pred[adjacent] = k;
	      } else {
		if (D[adjacent] > value[i] + D[k]) {
		  D[adjacent] = value[i] + D[k];
		  IndexArcPred[adjacent] = i;
		  //Pred[adjacent] = k + 1;
		  Pred[adjacent] = k;
		}
	      }
	    }
	  }
	  
	  /* Retire le noeud traite */
	  if (IndexMin != 0) {
	    IndexListNode[IndexMinPred] = IndexListNode[IndexMin];
	    if (IndexListNode[IndexMinPred] == n)
	      MaxNode = IndexMinPred + 1;
	  }
	} else {
	  j = n;
	}
      }
      /* Compute the gradient of commodity */
     
      for(unsigned int i = 0; i < nIndCom ; i++) {
	unsigned int k = _data->_dlistAdj[i + minIndCom];
	if (s != k) {
	  double demand = demandV[i + minIndCom];
	  grad[IndexArcPred[k] + _data->_numArcs * CurrentSourceNode] -= demand;
	  while(Pred[k] != s) {
	    unsigned int r = Pred[k];
	    grad[IndexArcPred[r] + _data->_numArcs * CurrentSourceNode] -= demand;
	    k = r;
	  }
	}
      }
      if(_aggreg == 2) {
	CurrentSourceNode += 1;
      }
    }
  }
  delete [] Pred;
  delete [] StateNode;
  delete [] D;
  delete [] IndexArcPred;
  delete [] ListNode;
  delete [] IndexListNode;

  return 0;
}

int 
MCFOracleFunction::eval(const AccpmVector &y, AccpmVector &functionValue, 
			AccpmGenMatrix &subGradients, AccpmGenMatrix *info)
{
  AccpmVector value(y);
  
  if (_data->_type == MCFData::LINEAR) {
    AccpmLAAddMult(value, 1, _data->_cost);
  }
  
  subGradients = 0;
  dijkstra(value.addr(), subGradients.addr());
   
  AccpmLAMatTransVecMult(subGradients, value, functionValue, 1, 0);
  
  if (info) {
    *info = AccpmGenMatrix(subGradients.size(1), 1);
    for (int i = 0; i < subGradients.size(1); ++i) {
      (*info)(i, 0) = i + 1;
    }
  }

  return _data->_stop;
}

const AccpmVector&
MCFOracleFunction::getCapacity() const
{
  return _data->_capacity;
}

int
MCFOracleFunction::updateAx(const AccpmVector &sol, double sumX) 
{
  _data->_sol = sol;
  _data->_sumX = sumX;
  return 0;
}

void 
MCFOracleFunction::setAggregationLevel(int level)
{
  if (level == 1 || level == 2) {
    _aggreg = level;
  } else {
    std::cerr << "Aggregation level " << level << " not supported. Defaulting to fully aggregated 1." << std::endl;
    _aggreg = 1;
  }
}
