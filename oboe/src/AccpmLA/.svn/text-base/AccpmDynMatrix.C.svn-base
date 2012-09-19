// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#include "AccpmDynMatrix.h"
#include <algorithm>
#include <iostream>

namespace Accpm 
{
static int POOLSIZE=100;

AccpmDynMatrix::AccpmDynMatrix(bool preallocateRows, bool preallocateCols, int m, int n) 
  : _preallocateRows(preallocateRows), _preallocateCols(preallocateCols),
    _currentRow(m), _currentCol(n), 
    _matrix(std::max(m, (preallocateRows * POOLSIZE)), std::max(n, (preallocateCols * POOLSIZE)))
  
{
}

AccpmDynMatrix::AccpmDynMatrix(const AccpmGenMatrix &rhs, bool preallocateRows, bool preallocateCols) 
  : _preallocateRows(preallocateRows), _preallocateCols(preallocateCols),
    _currentRow(rhs.size(0)), _currentCol(rhs.size(1)), 
    _matrix(rhs)
{
}

AccpmDynMatrix::AccpmDynMatrix(const AccpmDynMatrix &rhs)
  : _preallocateRows(rhs._preallocateRows), _preallocateCols(rhs._preallocateCols),
    _currentRow(rhs._currentRow), _currentCol(rhs._currentCol), _matrix(rhs._matrix)
{
}
  
AccpmDynMatrix& 
AccpmDynMatrix::operator=(double s)
{
  _matrix = s;
  return *this;
}

AccpmDynMatrix& 
AccpmDynMatrix::operator=(const AccpmDynMatrix &rhs)
{
  if (this == &rhs) {
    return *this;
  }
  _preallocateRows = rhs._preallocateRows;
  _preallocateCols = rhs._preallocateCols;
  _currentRow = rhs._currentRow;
  _currentCol = rhs._currentCol;
  _matrix = rhs._matrix;

  return *this;
}

AccpmDynMatrix::~AccpmDynMatrix()
{
}

AccpmDynMatrix&  
AccpmDynMatrix::resize(int m, int n)
{
  //std::cerr << "AccpmDynMat::Resizing"  << std::endl;
  RealMatrix tmp(m,n);
  if (_matrix.size(0) > 0 &&  _matrix.size(1) > 0) {
    if (m >= _matrix.size(0) && n >= _matrix.size(1)) {
      tmp(_matrix.index(0), _matrix.index(1)).inject(_matrix);
    } 
  }
  _matrix.ref(tmp);
  
  return *this;
}

int
AccpmDynMatrix::addColumn(const AccpmVector &v)
{
  if (_currentRow == 0) { // empty matrix
    _currentRow = v.size();
    if (_preallocateCols) {
      resize(_currentRow, _matrix.size(1) + POOLSIZE);
    } else {
      resize(_currentRow, _matrix.size(1) + 1);
    }
  }
  if (_currentCol == _matrix.size(1)) {
    if (_preallocateCols) {
      resize(_matrix.size(0), _matrix.size(1) + POOLSIZE);
    } else {
      resize(_matrix.size(0), _matrix.size(1) + 1);
    }
  }
  _matrix.assignColumn(_currentCol, v);
  ++_currentCol;
  return _currentCol-1;
 
}
   
int 
AccpmDynMatrix::addRow(const AccpmVector &v)
{
  if (_currentCol == 0) { // empty matrix
    _currentCol = v.size();
    if (_preallocateRows) {
      resize(_matrix.size(0) + POOLSIZE, _currentCol);
    } else {
      resize(_matrix.size(0) + 1, _currentCol);
    }
  }
  if (_currentRow == _matrix.size(0)) {
    if (_preallocateRows) {
      resize(_matrix.size(0) + POOLSIZE, _matrix.size(1));
    } else {
      resize(_matrix.size(0) + 1, _matrix.size(1));
    }
  }
  _matrix.assignRow(_currentRow, v);
  ++_currentRow;
  return _currentRow-1;

}

AccpmGenMatrix
AccpmDynMatrix::getM() const
{
  AccpmGenMatrix subMatrix;
  if (_currentRow > 0 && _currentCol > 0) {
    LaIndex rowIndex(0, _currentRow-1);
    LaIndex colIndex(0, _currentCol-1);
    subMatrix.resize(_currentRow, _currentCol);
    subMatrix.inject(_matrix(rowIndex, colIndex));
  }
  return subMatrix;
}

int
AccpmDynMatrix::size(int n) const
{
  if (n == 0) {
    return _currentRow;
  } 
  if (n == 1) {
    return _currentCol;
  } 
  return -1;
}    

int
AccpmDynMatrix::deleteColumn(int id)
{
  assert(0 <= id && id < size(1));
  AccpmGenMatrix tmp(_matrix.size(0), _matrix.size(1));
  LaIndex rowIndex(0, _currentRow-1);
 
  for (int i = 0; i < _currentCol; ++i) {
    LaIndex colIndex(i, i);
    if (i < id) {
      tmp(rowIndex, colIndex).inject(_matrix(rowIndex, colIndex));
    }
    if (i > id) {
      LaIndex tmpColIndex(i-1, i-1);
      tmp(rowIndex, tmpColIndex).inject(_matrix(rowIndex, colIndex)); 
    }
  }
  _matrix.ref(tmp);
  --_currentCol;
  
  return id;
}

}   
