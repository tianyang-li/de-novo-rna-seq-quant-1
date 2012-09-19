// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "Oracle.h"
#include "QpGenerator.h"
#include "Parameters.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <queue>

using namespace Accpm;
using namespace std;
using std::string;

class APSOracleFunction : public OracleFunction {

private:
  struct flight {
    flight(int ix=-1, int oc=-1,int row=-1)
      : ix_(ix), oc_(oc), row_(row) 
    { }
    
    int ix_;
    int oc_;
    int row_;
  };
  
  struct cont_t {
    int contribution_;
    flight vol_;
  };
  
  struct column {
    int cost_;
    std::vector<cont_t> rot_;
  };

  std::map<int,flight> _flights;
  std::vector<column> _matriceA;
  int *_demands;
  int _v_dot_b;

  static void print_flights(std::pair<int, flight> vol)
  {
    std::cout << "Flight_"<<vol.first<<" oc_"<<vol.second.oc_<<'\t';
  }

  void readAPCRotation(const char *fileName, bool explicitForm)
  {
    std::ifstream fin(fileName);
    int SZ=1024;
    char buf[SZ];
    
    fin.getline(buf,SZ);
    std::string str1(buf);
    std::string str2("FLIGHTS");
    assert(str1 == str2);
    int row_count=0;
    while(fin.getline(buf,SZ) && std::string(buf) != std::string("PAIRINGS")) {
      int flight_id, overcover_cost;
      std::istringstream str(buf);
      str >> flight_id >> overcover_cost;
      _flights.insert(std::pair<int,flight>(flight_id, flight(flight_id, overcover_cost, row_count)));
      ++row_count;
    }
    //std::for_each(_flights.begin(), _flights.end(), print_flights);
    std::cout<<std::endl<<std::endl;
      
    // store all pairings
    while(fin.getline(buf,SZ)) {
      std::istringstream str(buf);
      int ix, row_0, num_contributions;
      column col;
      str >> ix >> row_0 >> col.cost_ >> num_contributions;
      assert(row_0 == 0 && num_contributions > 0); //pairing well read...
      for (int i = 0; i < num_contributions; ++i) {
	int fl_ix, coeff;
	str >> fl_ix >> coeff;
	if (coeff != -1) { //omit deadheads  
	  flight vol(_flights[fl_ix]);
	  cont_t vol_t;
	  assert(coeff == 1);
	  vol_t.contribution_ = coeff;
	  vol_t.vol_ = vol;
	  col.rot_.push_back(vol_t);
          col.cost_ += !explicitForm ? vol.oc_ : 0;
	}
      }
      assert(col.rot_.size() > 0);
      _matriceA.push_back(col);
    }
  }

  void 
  addExplicitOvercoversVars()
  {
    std::cout << "Adding Explicit Overcovers" << std::endl;
    // add -I surplus variables
    int i=0;
    for(std::map<int,flight>::iterator iter = _flights.begin(); iter != _flights.end(); ++iter) {
      std::pair<int,flight> p = *iter;
      column col;
      col.cost_ = p.second.oc_;
      cont_t vol_t;
      vol_t.contribution_ = -1;
      vol_t.vol_ = p.second;
      col.rot_.push_back(vol_t);
      _matriceA.push_back(col);
      ++i;
    }
    std::cout<<"matriceA.size() : "<< _matriceA.size()<<std::endl;
  }

  double scalar_int( const double *x1, const int *x2, int n)
  {
    if( n <= 0 ) return 0;
    double sum=0;
    for(int k=0; k < n; ++k) sum += (x1[k])*(x2[k]);
    return sum;
  }
  
  struct AttractiveColumn {
    AttractiveColumn(const double & rc, const int & ix, const int *a, const int & dim)
      : rc_(rc), ix_(ix), a_(dim) 
    { 
      for (int i=0; i < dim; ++i)
	a_[i]=a[i];
    }
    
    void print() const {
      std::cout << "attractive colix="<<ix_<<"  rc="<<rc_<< std::endl;
    }
  
    double rc_;
    int ix_;
    std::vector<int> a_;
  };
  
  class compare_AttractiveColumn_pq {
  public:
    int operator () (const AttractiveColumn & lhs, const AttractiveColumn & rhs)
    {
      return lhs.rc_ > rhs.rc_; //ascending order for priority_queue
    };
  };

  void 
  compute_function_and_gradient_experiment_P_attractives(const double *y, const int & dim,
							 const int &P,  // choose P columns
							 const int &wP, // with weight wP 
							 double &My_f, double *My_g)
  { 
   
    int a[dim];

    // compute My_g = wP*Sum(most P attractive columns k){a[k]} - b, where a[k] is matrix column k.
    for (int i=0; i < dim; ++i)
      My_g[i]=-_demands[i];
    //memset(My_g, 0, sizeof(double)*dim);

    // compute My_f = <-b.y> + wP*Sum(most P attractive columns k){-RC[k]},  where RC[k] = c(k) - <y,a[k]>.
    My_f = - scalar_int(y, _demands, dim);
    double myf_all=My_f;
    
    //most attractive columns sorted in ascending order
    std::priority_queue<AttractiveColumn, std::vector<AttractiveColumn>, compare_AttractiveColumn_pq> attractives;
  
    int count_attractives=0;
    for (unsigned int k=0; k < _matriceA.size(); ++k) {
      memset(a, 0, sizeof(int)*dim);

      column col = _matriceA[k];
      int nnz = col.rot_.size();
      for (int i = 0; i < nnz; ++i) 
	a[col.rot_[i].vol_.row_] = col.rot_[i].contribution_;
      double reduced_cost = double(col.cost_) - scalar_int(y,a,dim);
      
      if (reduced_cost < 0) {
	attractives.push(AttractiveColumn(reduced_cost, k, a, dim));
	++count_attractives;
	myf_all += (-reduced_cost);
      }
      
#if 0
      double min_reduced_cost = FLT_MAX;
      int best_reduced_ix=-1;
      if (reduced_cost < min_reduced_cost) {
	min_reduced_cost = reduced_cost;
	best_reduced_ix = k;
      }
#endif
    }

    int cnt=0;
    while (attractives.size() > 0 && cnt < P) {
      const AttractiveColumn & p = attractives.top();
      //p.print();
      My_f += (wP*(-p.rc_));
      for(int i=0; i < dim; ++i) My_g[i] += (wP*p.a_[i]);
      attractives.pop();
      ++cnt;
    }
   
#if 0
    double min_dual = FLT_MAX;
    double max_dual = -FLT_MAX;
    for (int i=0; i < dim; ++i) {
      if (y[i] < min_dual) min_dual = y[i];
      if (y[i] > max_dual) max_dual = y[i];
    }
    printf("\n          min_dual=%lf  max_dual=%lf  min(rc)=%lf  col_ix=%d  attractives=%d\n", 
	   min_dual, max_dual, min_reduced_cost,best_reduced_ix,count_attractives);
    std::cout<<"     \t\t\t\t\t      "<<My_f<<"=L(y)   "<<"myf_all="<<myf_all<<std::endl;
#endif
  }

public:
  APSOracleFunction() : OracleFunction() {}
  APSOracleFunction(const char *fileName, bool explicitForm) : OracleFunction() {
    readAPCRotation(fileName, explicitForm);
    _demands = new int[_flights.size()];
    for (unsigned int i = 0; i < _flights.size(); ++i) {
      _demands[i] = 1;
    }
    _v_dot_b = 0;
    if (explicitForm) {
      addExplicitOvercoversVars();
    } else {
      for (std::map<int,flight>::iterator iter = _flights.begin(); 
	   iter != _flights.end(); ++iter) {
	assert((*iter).second.oc_ >= 0);
	_v_dot_b += (*iter).second.oc_;
      }
      std::cout<<"overcover_offset = <v,b=1> = "<< _v_dot_b << std::endl;
    }
  }
  virtual ~APSOracleFunction() { delete [] _demands; }
  
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info) {
    double fVal;
    double g[y.size()];
    int p = 100000;
    int wp = 1;
    compute_function_and_gradient_experiment_P_attractives(y.addr(), y.size(), p, wp, fVal, g);

    if (_v_dot_b > 0) { // incase of implicit formulation
      //std::cout<<"F + v_dot_b = " << fVal +_v_dot_b << std::endl;
      fVal += _v_dot_b;
    }

    functionValue = fVal;
    memcpy(subGradients.addr(), g, sizeof(double)*y.size());

    *info = 1;
    return 0;
  }
  
  static int 
  getNumberFlights(const char *fileName)
  {
    
    std::ifstream fin(fileName);
    int row_count = 0;
    if (fin.is_open()) {
      int SZ = 1024;
      char buf[SZ];
      fin.getline(buf,SZ);
      std::string str1(buf);
      std::string str2("FLIGHTS");
      assert(str1 == str2);
      
      while(fin.getline(buf,SZ) && std::string(buf) != std::string("PAIRINGS"))
	++row_count;
      fin.close();
      return row_count;
    } else {
      std::cout << "Error opening file: " << fileName << std::endl;
    }
    return row_count;
  }
};

int 
main(int argc, char *argv[])
{
  char *paramFile = "param.txt";
  char *rotationFile = "APC_rot.reservoir";
  bool explicitForm = true;
  if (argc != 3 && argc != 4) {
    std::cout << "Error starting " << argv[0] << " with " << argc  << " arguments" <<std::endl;
    std::cout << "Usage: " << argv[0] << " <Parameter File> <Rotation File> <explicit formualtion>"
	<< std::endl;
    exit(0);
  }
  paramFile = argv[1];
  rotationFile = argv[2];
  if (argc == 4) {
    explicitForm = atoi(argv[3]);
  }
  
  Accpm::Parameters param(paramFile);
  int n = APSOracleFunction::getNumberFlights(rotationFile);
  param.setIntParameter("NumVariables",n);

  vector<double> start(n, 200);
  param.setStartingPoint(start);
  //vector<double> b(n, -1);
  //param.setB(b);
  if (!explicitForm) {
    std::cout << "Using implicit formulation" << std::endl;
    vector<double> varLB(n, 0);
    param.setVariableLB(varLB);
  }

  APSOracleFunction f1(rotationFile, explicitForm);
  Accpm::Oracle oracle(&f1);
  
  QpGenerator qpGen;
  qpGen.init(&param, &oracle);
 
  while (!qpGen.run()) {

  }
  qpGen.output(cout);
  qpGen.terminate();
}

