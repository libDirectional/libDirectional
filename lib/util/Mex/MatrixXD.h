
/* >> This file is part of the Nonlinear Estimation Toolbox
 *
 *    For more information, see https://bitbucket.org/NonlinearEstimation/toolbox
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

#ifndef _MEX_MATRIX_XD_H_
#define _MEX_MATRIX_XD_H_

#include "Traits.h"
#include "Utils.h"
#include <iostream>

#define isEigenAligned                                                                  \
    ((r != Eigen::Dynamic && r * sizeof(Scalar) % 16 == 0) ||                           \
     (c != Eigen::Dynamic && c * sizeof(Scalar) % 16 == 0) ||                           \
     (r != Eigen::Dynamic && c != Eigen::Dynamic && r * c * sizeof(Scalar) % 16 == 0)	\
     ? Eigen::Aligned : Eigen::Unaligned)                                               \

namespace Mex {

// iterator for multidimensional slices.
// see struct SliceRange
class SliceIterator : public std::input_iterator_tag {
    public:
        SliceIterator(const SliceIterator& it) :
            pos(it.pos),
            dims(it.dims) { }
        
        explicit SliceIterator(const Eigen::VectorXi& d) {
            if (d.size() == 0) {
                pos  = Eigen::VectorXi::Zero(1, 1);
                dims = Eigen::VectorXi::Ones(1, 1);
            } else {
                pos  = Eigen::VectorXi::Zero(d.size());
                dims = d;
            }
        }
        
        virtual ~SliceIterator() { }
        
        static SliceIterator end(const Eigen::VectorXi& dims) {
            SliceIterator it(dims);
            
            it.pos(it.dims.size() - 1) = it.dims(it.dims.size() - 1);
            
            return it;
        }
        
        bool operator==(const SliceIterator& it) const {
            return pos == it.pos;
        }
        
        bool operator!=(const SliceIterator& it) const {
            return pos != it.pos;
        }
        
        const Eigen::VectorXi& operator*() const {
            return pos;
        }
        
        SliceIterator& operator++() {
            int i;
            
            for (i = 0; i < pos.size() - 1; i++) {
                pos(i)++;
                
                if (pos(i) == dims(i)) {
                    pos(i) = 0;
                } else {
                    return *this;
                }
            }
            
            pos(i)++;
            
            return *this;
        }
        
        SliceIterator operator++(int) {
            SliceIterator temp = *this;
            
            ++*this;
            
            return temp;
        }
        
    private:
        Eigen::VectorXi pos;
        Eigen::VectorXi dims;
        
};

// iterator for multidimensional slices. example use:
// ConstMatrixXDX mat{3, 4, 5, 6};
// for (const Eigen::VectorXi &index : SliceRange(mat.sliceSize()) {
//      Matrix slice = mat.slice(index);
//      ... // slice has size (3, 4); index ranges from (0, 0) to (5, 6)
// }
class SliceRange {
    public:
        template<class Mat>
        SliceRange(const Mat& mat) : dims(mat.slices()) { }
        
        SliceRange(const Eigen::VectorXi& d) : dims(d) { }
        
        SliceIterator begin() const {
            return SliceIterator(dims);
        }
        
        SliceIterator end() const {
            return SliceIterator::end(dims);
        }
        
    private:
        const Eigen::VectorXi dims;
        
};

template<typename Scalar, int r, int c>
class ConstMatrixXD {
    public:
        template<bool ForceAlign> 
        using ConstSliceBase = Eigen::Map<const Eigen::Matrix<Scalar, r, c>, 
                                          ForceAlign ? Eigen::Aligned : isEigenAligned>;
        
        typedef ConstSliceBase<false> ConstSlice;
        
    public:
        ConstMatrixXD(const mxArray* array) :
            mxData(array),
            matData(Utils::checkArrayType<Scalar>(array)),
            dimensions(Utils::getDimensions(array)) {
          	Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            setArrayData();
        }
        
        int rows() const {
            return dimensions(0);
        }
        
        int cols() const {
            return dimensions(1);
        }
        
        int dim(unsigned int i) const {
            return i < dimensions.size() ? dimensions(i) : 1;
        }
        
        const Eigen::VectorXi& dims() const {
            return dimensions;
        }
        
        const Eigen::VectorXi& slices() const {
            return sliceDims;
        }
        
        int size() const {
            return numElements;
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(const Eigen::VectorXi& reqSlice) const {
            mxAssert(Utils::isValidSlice(reqSlice, sliceMins, sliceMaxs),
                     "Requested slice is out of range");
            
            const int offset = reqSlice.head(sliceOffsets.size()).dot(sliceOffsets);
            
        	const Scalar* sliceData = &matData[offset];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(int reqSlice) const {
            const int element = reqSlice * sliceSize;
            
            mxAssert(element < numElements && element >= 0,
                     "Requested slice is out of range");
            
        	const Scalar* sliceData = &matData[element];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> sliceExpanded(const Eigen::VectorXi& reqSlice) const {
            Eigen::VectorXi expSlice = reqSlice.head(sliceDims.size());
            
            for (int i = 0; i < sliceDims.size(); i++) {
                if (sliceDims(i) == 1) {
                    expSlice(i) = 0;
                }
            }
            
            return ConstSliceBase<ForceAlign>(expSlice);
        }
        
        void print() const {
            for (const Eigen::VectorXi& s : SliceRange(*this)) {
                std::cout << "this(:,:," << s.transpose() << ") =\n"
                          << slice(s) << std::endl << std::endl;
            }
        }
        
        const Scalar* data() const {
            return matData;
        }
        
        operator const mxArray*() const {
            return mxData;
        }
        
    private:
        void setArrayData() {
            sliceDims = dimensions.tail(dimensions.size() - 2);
            
            sliceOffsets.resize(sliceDims.size());
            
            sliceSize = rows() * cols();
            
            int start = sliceSize;
            
            for (int i = 0; i < sliceOffsets.size(); i++) {
                sliceOffsets(i) = start;
                start *= dimensions(i + 2);
            }
            
            numElements = start;
            
            sliceMins = Eigen::VectorXi::Zero(sliceDims.size());
            sliceMaxs = sliceDims.array() - 1;
        }
        
        /* Non-copyable */
        ConstMatrixXD(const ConstMatrixXD& mat) = delete;
        ConstMatrixXD& operator=(const ConstMatrixXD& mat) = delete;
        
    private:
        const mxArray*     	const mxData;
        const Scalar*    	const matData;
        Eigen::VectorXi     dimensions;
        Eigen::VectorXi   	sliceDims;
        Eigen::VectorXi  	sliceOffsets;
        Eigen::VectorXi   	sliceMins;
        Eigen::VectorXi    	sliceMaxs;
        int                 numElements;
        int                 sliceSize;
        
};

template<bool hasOwnership, typename Scalar, int r, int c>
class MatrixBaseXD {
    public:
        template< bool ForceAlign> 
        using ConstSliceBase = Eigen::Map<const Eigen::Matrix<Scalar, r, c>, 
                                          ForceAlign ? Eigen::Aligned : isEigenAligned>;
        
     	typedef ConstSliceBase<false> ConstSlice;
        
        template<bool ForceAlign> 
        using SliceBase = Eigen::Map<Eigen::Matrix<Scalar, r, c>, 
                                     ForceAlign ? Eigen::Aligned : isEigenAligned>;
        
     	typedef SliceBase<false> Slice;
        
    public:
        // Usage: MatrixBaseXD mat{3, 4, 5, 6, 7};
        MatrixBaseXD(std::initializer_list<int> reqDims) :
            dimensions(createDims(reqDims)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            createData();
        }
        
        MatrixBaseXD(const Eigen::VectorXi& reqDims) :
            dimensions(createDims(reqDims)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            createData();
        }
        
        MatrixBaseXD(int numRows, int numCols,
                     const Eigen::VectorXi& reqSliceDims) :
            dimensions(createDims(numRows, numCols, reqSliceDims)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            createData();
     	}
        
        MatrixBaseXD(mxArray* array) :
            mxData(array),
            matData(Utils::checkArrayType<Scalar>(array)),
            dimensions(Utils::getDimensions(array)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            setArrayData();
        }
        
        virtual ~MatrixBaseXD() {
            if (hasOwnership) {
                mxDestroyArray(mxData);
            }
        }
        
        int rows() const {
            return dimensions(0);
        }
        
        int cols() const {
            return dimensions(1);
        }
        
        int dim(unsigned int i) const {
            return i < dimensions.size() ? dimensions(i) : 1;
        }
        
        const Eigen::VectorXi& dims() const {
            return dimensions;
        }
        
        const Eigen::VectorXi& slices() const {
            return sliceDims;
        }
        
        int size() const {
            return numElements;
        }
        
        template<bool ForceAlign = false>
        SliceBase<ForceAlign> slice(const Eigen::VectorXi& reqSlice) {
            mxAssert(Utils::isValidSlice(reqSlice, sliceMins, sliceMaxs),
                     "Requested slice is out of range");
            
            const int offset = reqSlice.head(sliceOffsets.size()).dot(sliceOffsets);
            
        	Scalar* sliceData = &matData[offset];
            
        	return SliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        SliceBase<ForceAlign> slice(int reqSlice) {
            const int element = reqSlice * sliceSize;
            
            mxAssert(element < numElements && element >= 0,
                     "Requested slice is out of range");
            
            Scalar* sliceData = &matData[element];
            
        	return SliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(const Eigen::VectorXi& reqSlice) const {
            mxAssert(Utils::isValidSlice(reqSlice, sliceMins, sliceMaxs),
                     "Requested slice is out of range");
            
            const int offset = reqSlice.head(sliceOffsets.size()).dot(sliceOffsets);
            
        	const Scalar* sliceData = &matData[offset];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(int reqSlice) const {
            const int element = reqSlice * sliceSize;
            
            mxAssert(element < numElements && element >= 0,
                     "Requested slice is out of range");
            
        	const Scalar* sliceData = &matData[element];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        SliceBase<ForceAlign> sliceExpanded(const Eigen::VectorXi& reqSlice) {
            Eigen::VectorXi expSlice = reqSlice.head(sliceDims.size());
            
            for (int i = 0; i < sliceDims.size(); i++) {
                if (sliceDims(i) == 1) {
                    expSlice(i) = 0;
                }
            }
            
            return SliceBase<ForceAlign>(expSlice);
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> sliceExpanded(const Eigen::VectorXi& reqSlice) const {
            Eigen::VectorXi expSlice = reqSlice.head(sliceDims.size());
            
            for (int i = 0; i < sliceDims.size(); i++) {
                if (sliceDims(i) == 1) {
                    expSlice(i) = 0;
                }
            }
            
            return ConstSliceBase<ForceAlign>(expSlice);
        }
        
        void print() const {
            for (const Eigen::VectorXi& s : SliceRange(slices())) {
                std::cout << "this(:,:," << s.transpose() << ") =\n"
                          << slice(s) << std::endl << std::endl;
            }
        }
        
        Scalar* data() {
            return matData;
        }
        
        const Scalar* data() const {
            return matData;
        }
        
        operator mxArray*() {
            return mxData;
        }
        
       	operator const mxArray*() const {
            return mxData;
        }
        
    private:
        Eigen::VectorXi createDims(const std::initializer_list<int>& reqDims) {
        	if (reqDims.size() < 1) {
                throw std::invalid_argument("At least one dimension must be provided.");
            } else if (reqDims.size() == 1) {
                Eigen::VectorXi vec(2);
                
                std::copy(reqDims.begin(), reqDims.end(), vec.data());
                
                vec(1) = 1;
                
                return vec;
            } else {
                Eigen::VectorXi vec(reqDims.size());
                
                std::copy(reqDims.begin(), reqDims.end(), vec.data());
                
                return vec;
            }
        }
        
        Eigen::VectorXi createDims(const Eigen::VectorXi& reqDims) {
        	if (reqDims.size() < 1) {
                throw std::invalid_argument("At least one dimension must be provided.");
            } else if (reqDims.size() == 1) {
                Eigen::VectorXi vec(2);
                
                vec(0) = reqDims(0);                
                vec(1) = 1;
                
                return vec;
            } else {
                return reqDims;
            }
        }
        
     	Eigen::VectorXi createDims(int rows, int cols, 
                                   const Eigen::VectorXi& reqSliceDims) {
            Eigen::VectorXi vec(2 + reqSliceDims.size());
            
            vec(0) = rows;
            vec(1) = cols;
            
            vec.tail(reqSliceDims.size()) = reqSliceDims;
            
            return vec;
        }
        
        void createData() {
            mxData  = mxCreateNumericArray(dimensions.size(), dimensions.data(), 
                                          Traits<Scalar>::MxClassID,  mxREAL);
            matData = (Scalar*)mxGetData(mxData);
            
            setArrayData();
        }
        
        void setArrayData() {
            sliceDims = dimensions.tail(dimensions.size() - 2);
            
            sliceOffsets.resize(sliceDims.size());
            
            sliceSize = rows() * cols();
            
            int start = sliceSize;
            
            for (int i = 0; i < sliceOffsets.size(); i++) {
                sliceOffsets(i) = start;
                start *= dimensions(2 + i);
            }
            
            numElements = start;
            
            sliceMins = Eigen::VectorXi::Zero(sliceDims.size());
            sliceMaxs = sliceDims.array() - 1;
        }
        
        /* Non-copyable */
        MatrixBaseXD(const MatrixBaseXD& mat) = delete;
        template<bool o, typename s, int rm, int cm>
        MatrixBaseXD(const MatrixBaseXD<o, s, rm, cm>& mat) = delete;
        MatrixBaseXD& operator=(const MatrixBaseXD& mat) = delete;
        
    private:
        mxArray*          	mxData;
        Scalar*          	matData;
      	Eigen::VectorXi     dimensions;
        Eigen::VectorXi   	sliceDims;
        Eigen::VectorXi   	sliceOffsets;
        Eigen::VectorXi    	sliceMins;
        Eigen::VectorXi     sliceMaxs;
        int                 numElements;
        int                 sliceSize;
        
};

template<typename Scalar = double, int r = Eigen::Dynamic, int c = Eigen::Dynamic>
using MatrixXD = MatrixBaseXD<true, Scalar, r, c>;

template<typename Scalar = double, int r = Eigen::Dynamic, int c = Eigen::Dynamic>
using OutputMatrixXD = MatrixBaseXD<false, Scalar, r, c>;


/* ConstMatrix typedefs */
template<typename Scalar = double>
using ConstMatrixXDX = ConstMatrixXD<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar = double>
using ConstVectorXDX = ConstMatrixXD<Scalar, Eigen::Dynamic, 1>;

template<typename Scalar = double>
using ConstRowVectorXDX = ConstMatrixXD<Scalar, 1, Eigen::Dynamic>;


/* Matrix typedefs */
template<typename Scalar = double>
using MatrixXDX = MatrixXD<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar = double>
using VectorXDX = MatrixXD<Scalar, Eigen::Dynamic, 1>;

template<typename Scalar = double>
using RowVectorXDX = MatrixXD<Scalar, 1, Eigen::Dynamic>;


/* OutputMatrix typedefs */
template<typename Scalar = double>
using OutputMatrixXDX = OutputMatrixXD<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar = double>
using OutputVectorXDX = OutputMatrixXD<Scalar, Eigen::Dynamic, 1>;

template<typename Scalar = double>
using OutputRowVectorXDX = OutputMatrixXD<Scalar, 1, Eigen::Dynamic>;


/* Double typedefs */
typedef ConstMatrixXDX<>        ConstMatrixXDXd;
typedef ConstVectorXDX<>        ConstVectorXDXd;
typedef ConstRowVectorXDX<>     ConstRowVectorXDXd;

typedef MatrixXDX<>             MatrixXDXd;
typedef VectorXDX<>             VectorXDXd;
typedef RowVectorXDX<>          RowVectorXDXd;

typedef OutputMatrixXDX<>    	OutputMatrixXDXd;
typedef OutputVectorXDX<>       OutputVectorXDXd;
typedef OutputRowVectorXDX<>    OutputRowVectorXDXd;

}   // namespace Mex

#endif
