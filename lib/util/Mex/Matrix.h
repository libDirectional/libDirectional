
/* >> This file is part of the Nonlinear Estimation Toolbox
 *
 *    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
 *
 *    Copyright (C) 2015  Jannik Steinbring <jannik.steinbring@kit.edu>
 *                        Antonio Zea <antonio.zea@kit.edu>  
 *
 *                        Institute for Anthropomatics and Robotics
 *                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
 *                        Karlsruhe Institute of Technology (KIT), Germany
 *
 *                        http://isas.uka.de
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MEX_MATRIX_H_
#define _MEX_MATRIX_H_

#include "Traits.h"
#include "Utils.h"
#include <iostream>

namespace Mex {

template<typename Scalar, int r, int c>
class ConstMatrix : public Eigen::Map<const Eigen::Matrix<Scalar, r, c>, Eigen::Aligned> { 
    private:
        typedef Eigen::Map<const Eigen::Matrix<Scalar, r, c>, Eigen::Aligned>  BaseClass;
        
    public:
        typedef typename Eigen::Matrix<Scalar, r, c>::Index  Index;
        
    public:
        ConstMatrix(const mxArray* array) :
            BaseClass(Utils::checkArray<Scalar, r, c>(array),
                      mxGetM(array),
                      mxGetN(array)),
            mxData(array) { }
      	
        void print() const {
            std::cout << *this << std::endl << std::endl;
        }
        
        operator const mxArray*() const {
            return mxData;
        }
        
    private:
        /* Non-copyable */
        ConstMatrix(const ConstMatrix& mat);
        ConstMatrix& operator=(const ConstMatrix& mat);
        
    private:
        const mxArray* const mxData;
        
};

template<bool hasOwnership, typename Scalar, int r, int c>
class MatrixBase : public Eigen::Map<Eigen::Matrix<Scalar, r, c>, Eigen::Aligned> {
    private:
        typedef Eigen::Map<Eigen::Matrix<Scalar, r, c>, Eigen::Aligned>  BaseClass;
        
    public:
        typedef typename Eigen::Matrix<Scalar, r, c>::Index  Index;
        
    public:
        MatrixBase() :
            BaseClass(createMxArray(),
                      r,
                      c) { }
      	
        MatrixBase(Index dim) :
            BaseClass(createVectorMxArray(dim),
                      r == Eigen::Dynamic ? dim : r,
                      c == Eigen::Dynamic ? dim : c) { }
    	
        MatrixBase(Index numRows,
                   Index numCols) :
            BaseClass(createMxArray(numRows, numCols),
                      numRows,
                      numCols) { }
        
        MatrixBase(mxArray* array) :
            BaseClass(Utils::checkArray<Scalar, r, c>(array),
                      mxGetM(array),
                      mxGetN(array)),
            mxData(array) { }
      	
        MatrixBase(const MatrixBase& mat) :
            BaseClass(createMxArray(mat.rows(), mat.cols()),
                      mat.rows(),
                      mat.cols()) {
            BaseClass::operator=(mat);
        }
     	
        template<typename Derived>
        MatrixBase(const Eigen::MatrixBase<Derived>& mat) :
            BaseClass(createMxArray(mat.rows(), mat.cols()),
                      mat.rows(),
                      mat.cols()) {
            BaseClass::operator=(mat); 
        }
        
        ~MatrixBase() {
            if (hasOwnership) {
                mxDestroyArray(mxData);
            }
        }
        
        MatrixBase& operator=(const MatrixBase& mat) {
            mxAssert(mat.rows() == this->rows(), "Different number of rows");
            mxAssert(mat.cols() == this->cols(), "Different number of cols");
            
            BaseClass::operator=(mat);
            
            return *this;
        }
        
        template<typename Derived>
        MatrixBase& operator=(const Eigen::MatrixBase<Derived>& mat) {
            mxAssert(mat.rows() == this->rows(), "Different number of rows");
            mxAssert(mat.cols() == this->cols(), "Different number of cols");
            
            BaseClass::operator=(mat);
            
            return *this;
        }
        
        void print() const {
            std::cout << *this << std::endl << std::endl;
        }
        
        operator mxArray*() {
            return mxData;
        }
        
        operator const mxArray*() const {
            return mxData;
        }
        
    private:
        Scalar* createMxArray() {
            static_assert(r != Eigen::Dynamic && c != Eigen::Dynamic,
                          "Both template parameters must be non-dynamic.");
            
            mxData = mxCreateNumericMatrix(r, c, Traits<Scalar>::MxClassID, mxREAL);
            
            return (Scalar*)mxGetData(mxData);
        }
        
        Scalar* createMxArray(Index rows,
                              Index cols) {
            Utils::checkRows<r>(rows);
            Utils::checkCols<c>(cols);
            
            mxData = mxCreateNumericMatrix(rows, cols, Traits<Scalar>::MxClassID, mxREAL);
            
            return (Scalar*)mxGetData(mxData);
        }
        
        Scalar* createVectorMxArray(Index dim) {
            static_assert(r == 1 || c == 1,
                          "You called a vector method an a matrix.");
            
            if (r == 1) {
                Utils::checkCols<c>(dim);
                
                mxData = mxCreateNumericMatrix(1, dim, Traits<Scalar>::MxClassID, mxREAL);
            } else {
                Utils::checkRows<r>(dim);
                
                mxData = mxCreateNumericMatrix(dim, 1, Traits<Scalar>::MxClassID, mxREAL);
            }
            
            return (Scalar*)mxGetData(mxData);
        }
        
    private:
        mxArray* mxData;
        
};


/* ConstMatrix typedefs */
template<typename Scalar = double>
using ConstMatrixX = ConstMatrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, int r>
using ConstVector = ConstMatrix<Scalar, r, 1>;

template<typename Scalar = double>
using ConstVectorX = ConstVector<Scalar, Eigen::Dynamic>;

template<typename Scalar, int c>
using ConstRowVector = ConstMatrix<Scalar, 1, c>;

template<typename Scalar = double>
using ConstRowVectorX = ConstRowVector<Scalar, Eigen::Dynamic>;


/* Matrix typedefs */
template<typename Scalar, int r, int c>
using Matrix = MatrixBase<true, Scalar, r, c>;

template<typename Scalar = double>
using MatrixX = Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, int r>
using Vector = Matrix<Scalar, r, 1>;

template<typename Scalar = double>
using VectorX = Vector<Scalar, Eigen::Dynamic>;

template<typename Scalar, int c>
using RowVector = Matrix<Scalar, 1, c>;

template<typename Scalar = double>
using RowVectorX = RowVector<Scalar, Eigen::Dynamic>;


/* OutputMatrix typedefs */
template<typename Scalar, int r, int c>
using OutputMatrix = MatrixBase<false, Scalar, r, c>;

template<typename Scalar = double>
using OutputMatrixX = OutputMatrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, int r>
using OutputVector = OutputMatrix<Scalar, r, 1>;

template<typename Scalar = double>
using OutputVectorX = OutputVector<Scalar, Eigen::Dynamic>;

template<typename Scalar, int c>
using OutputRowVector = OutputMatrix<Scalar, 1, c>;

template<typename Scalar = double>
using OutputRowVectorX = OutputRowVector<Scalar, Eigen::Dynamic>;


/* Double typedefs */
typedef ConstMatrixX<>      ConstMatrixXd;
typedef MatrixX<>           MatrixXd;
typedef OutputMatrixX<>     OutputMatrixXd;

typedef ConstVectorX<>      ConstVectorXd;
typedef VectorX<>           VectorXd;
typedef OutputVectorX<> 	OutputVectorXd;

typedef ConstRowVectorX<>   ConstRowVectorXd;
typedef RowVectorX<>        RowVectorXd;
typedef OutputRowVectorX<>	OutputRowVectorXd;

}   // namespace Mex

#endif
