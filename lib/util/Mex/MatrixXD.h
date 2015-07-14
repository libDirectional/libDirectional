
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

#ifndef _MEX_MATRIX_XD_H_
#define _MEX_MATRIX_XD_H_

#include "Traits.h"
#include "Utils.h"
#include <iostream>

#define isEigenAligned                                                                      \
    ((r != Eigen::Dynamic && r * (int)sizeof(Scalar) % 16 == 0) ||                          \
     (c != Eigen::Dynamic && c * (int)sizeof(Scalar) % 16 == 0) ||                          \
     (r != Eigen::Dynamic && c != Eigen::Dynamic && r * c * (int)sizeof(Scalar) % 16 == 0)  \
     ? Eigen::Aligned : Eigen::Unaligned)                                                   \

namespace Mex {

// iterator for multidimensional slices.
// see struct SliceRange
class SliceIterator : public std::input_iterator_tag {
    public:
        SliceIterator(const SliceIterator& it) :
            pos(it.pos),
            dims(it.dims) { }
        
        explicit SliceIterator(const Dimensions& d) {
            if (d.size() == 0) {
                pos  = Dimensions::Zero(1, 1);
                dims = Dimensions::Ones(1, 1);
            } else {
                pos  = Dimensions::Zero(d.size());
                dims = d;
            }
        }
        
        virtual ~SliceIterator() { }
        
        static SliceIterator end(const Dimensions& dims) {
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
        
        const Dimensions& operator*() const {
            return pos;
        }
        
        SliceIterator& operator++() {
            int64_t i;
            
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
        Dimensions pos;
        Dimensions dims;
        
};

// iterator for multidimensional slices. example use:
// Mex::ConstMatrixXDX mat{3, 4, 5, 6};
// for (const Mex::Dimensions &index : Mex::SliceRange(mat.sliceSize()) {
//      Matrix slice = mat.slice(index);
//      ... // slice has size (3, 4); index ranges from (0, 0) to (5, 6)
// }
class SliceRange {
    public:
        template<class Mat>
        SliceRange(const Mat& mat) : dims(mat.slices()) { }
        
        SliceRange(const Dimensions& d) : dims(d) { }
        
        SliceIterator begin() const {
            return SliceIterator(dims);
        }
        
        SliceIterator end() const {
            return SliceIterator::end(dims);
        }
        
    private:
        const Dimensions dims;
        
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
            matData(Utils::checkArray<Scalar, r, c, false>(array)),
            dimensions(Utils::getDimensions(array)) {
            setArrayData();
        }
        
        int64_t rows() const {
            return dimensions(0);
        }
        
        int64_t cols() const {
            return dimensions(1);
        }
        
        int64_t dim(int64_t i) const {
            return i < dimensions.size() ? dimensions(i) : 1;
        }
        
        const Dimensions& dims() const {
            return dimensions;
        }
        
        const Dimensions& slices() const {
            return sliceDims;
        }
        
        int64_t size() const {
            return numElements;
        }
        
        const Scalar& operator()(int64_t index) const {
            return matData[index];
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(const Dimensions& reqSlice) const {
            mxAssert(Utils::isValidSlice(reqSlice, sliceMins, sliceMaxs),
                     "Requested slice is out of range");
            
            const int64_t offset = reqSlice.head(sliceOffsets.size()).dot(sliceOffsets);
            
        	const Scalar* sliceData = &matData[offset];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(int64_t reqSlice) const {
            const int64_t element = reqSlice * sliceSize;
            
            mxAssert(element < numElements && element >= 0,
                     "Requested slice is out of range");
            
        	const Scalar* sliceData = &matData[element];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> sliceExpanded(const Dimensions& reqSlice) const {
            Dimensions expSlice = reqSlice.head(sliceDims.size());
            
            for (int64_t i = 0; i < sliceDims.size(); i++) {
                if (sliceDims(i) == 1) {
                    expSlice(i) = 0;
                }
            }
            
            return ConstSliceBase<ForceAlign>(expSlice);
        }
        
        void print() const {
            for (const Dimensions& s : SliceRange(*this)) {
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
            
            int64_t start = sliceSize;
            
            for (int64_t i = 0; i < sliceOffsets.size(); i++) {
                sliceOffsets(i) = start;
                start *= dimensions(i + 2);
            }
            
            numElements = start;
            
            sliceMins = Dimensions::Zero(sliceDims.size());
            sliceMaxs = sliceDims.array() - 1;
        }
        
        /* Non-copyable */
        ConstMatrixXD(const ConstMatrixXD& mat) = delete;
        ConstMatrixXD& operator=(const ConstMatrixXD& mat) = delete;
        
    private:
        const mxArray*      const mxData;
        const Scalar*       const matData;
        Dimensions          dimensions;
        Dimensions          sliceDims;
        Dimensions          sliceOffsets;
        Dimensions          sliceMins;
        Dimensions          sliceMaxs;
        int64_t           	numElements;
        int64_t           	sliceSize;
        
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
        MatrixBaseXD(std::initializer_list<int64_t> reqDims) :
            dimensions(createDims(reqDims)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            createData();
        }
        
        MatrixBaseXD(const Dimensions& reqDims) :
            dimensions(createDims(reqDims)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            createData();
        }
        
        MatrixBaseXD(int64_t numRows,
                     int64_t numCols,
                     const Dimensions& reqSliceDims) :
            dimensions(createDims(numRows, numCols, reqSliceDims)) {
            Utils::checkRows<r>(rows());
           	Utils::checkCols<c>(cols());
            
            createData();
     	}
        
        MatrixBaseXD(mxArray* array) :
            mxData(array),
            matData(Utils::checkArray<Scalar, r, c, false>(array)),
            dimensions(Utils::getDimensions(array)) {
            setArrayData();
        }
        
        virtual ~MatrixBaseXD() {
            if (hasOwnership) {
                mxDestroyArray(mxData);
            }
        }
        
        int64_t rows() const {
            return dimensions(0);
        }
        
        int64_t cols() const {
            return dimensions(1);
        }
        
        int64_t dim(int64_t i) const {
            return i < dimensions.size() ? dimensions(i) : 1;
        }
        
        const Dimensions& dims() const {
            return dimensions;
        }
        
        const Dimensions& slices() const {
            return sliceDims;
        }
        
        int64_t size() const {
            return numElements;
        }
        
        void setOnes() {
            for (int64_t i = 0; i < size(); ++i) {
                matData[i] = 1;
            }
        }
        
        void setZero() {
            std::memset(matData, 0, sizeof(Scalar) * size());
        }
        
        Scalar& operator()(int64_t index) {
            return matData[index];
        }
        
        const Scalar& operator()(int64_t index) const {
            return matData[index];
        }
        
        template<bool ForceAlign = false>
        SliceBase<ForceAlign> slice(const Dimensions& reqSlice) {
            mxAssert(Utils::isValidSlice(reqSlice, sliceMins, sliceMaxs),
                     "Requested slice is out of range");
            
            const int64_t offset = reqSlice.head(sliceOffsets.size()).dot(sliceOffsets);
            
        	Scalar* sliceData = &matData[offset];
            
        	return SliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        SliceBase<ForceAlign> slice(int64_t reqSlice) {
            const int64_t element = reqSlice * sliceSize;
            
            mxAssert(element < numElements && element >= 0,
                     "Requested slice is out of range");
            
            Scalar* sliceData = &matData[element];
            
        	return SliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(const Dimensions& reqSlice) const {
            mxAssert(Utils::isValidSlice(reqSlice, sliceMins, sliceMaxs),
                     "Requested slice is out of range");
            
            const int64_t offset = reqSlice.head(sliceOffsets.size()).dot(sliceOffsets);
            
        	const Scalar* sliceData = &matData[offset];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> slice(int64_t reqSlice) const {
            const int64_t element = reqSlice * sliceSize;
            
            mxAssert(element < numElements && element >= 0,
                     "Requested slice is out of range");
            
        	const Scalar* sliceData = &matData[element];
            
        	return ConstSliceBase<ForceAlign>(sliceData, rows(), cols());
        }
        
        template<bool ForceAlign = false>
        SliceBase<ForceAlign> sliceExpanded(const Dimensions& reqSlice) {
            Dimensions expSlice = reqSlice.head(sliceDims.size());
            
            for (int64_t i = 0; i < sliceDims.size(); i++) {
                if (sliceDims(i) == 1) {
                    expSlice(i) = 0;
                }
            }
            
            return SliceBase<ForceAlign>(expSlice);
        }
        
        template<bool ForceAlign = false>
        ConstSliceBase<ForceAlign> sliceExpanded(const Dimensions& reqSlice) const {
            Dimensions expSlice = reqSlice.head(sliceDims.size());
            
            for (int64_t i = 0; i < sliceDims.size(); i++) {
                if (sliceDims(i) == 1) {
                    expSlice(i) = 0;
                }
            }
            
            return ConstSliceBase<ForceAlign>(expSlice);
        }
        
        void print() const {
            for (const Dimensions& s : SliceRange(slices())) {
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
        Dimensions createDims(const std::initializer_list<int64_t>& reqDims) {
            if (reqDims.size() < 2) {
                throw std::invalid_argument("At least two dimensions must be provided.");
            } else {
                Dimensions vec(reqDims.size());
                
                std::copy(reqDims.begin(), reqDims.end(), vec.data());
                
                return vec;
            }
        }
        
        Dimensions createDims(const Dimensions& reqDims) {
            if (reqDims.size() < 2) {
                throw std::invalid_argument("At least two dimensions must be provided.");
            } else {
                return reqDims;
            }
        }
        
        Dimensions createDims(int64_t rows,
                              int64_t cols,
                              const Dimensions& reqSliceDims) {
            Dimensions vec(2 + reqSliceDims.size());
            
            vec(0) = rows;
            vec(1) = cols;
            
            vec.tail(reqSliceDims.size()) = reqSliceDims;
            
            return vec;
        }
        
        void createData() {
            const mwSize numDims = dimensions.size();
            mwSize* dims = new mwSize[numDims];
            
            for (int64_t i = 0; i < dimensions.size(); ++i) {
                dims[i] = dimensions[i];
            }
            
            mxData = mxCreateNumericArray(numDims, dims, Traits<Scalar>::MxClassID,  mxREAL);
            
            delete[] dims;
            
            matData = (Scalar*)mxGetData(mxData);
            
            setArrayData();
        }
        
        void setArrayData() {
            sliceDims = dimensions.tail(dimensions.size() - 2);
            
            sliceOffsets.resize(sliceDims.size());
            
            sliceSize = rows() * cols();
            
            int64_t start = sliceSize;
            
            for (int64_t i = 0; i < sliceOffsets.size(); i++) {
                sliceOffsets(i) = start;
                start *= dimensions(2 + i);
            }
            
            numElements = start;
            
            sliceMins = Dimensions::Zero(sliceDims.size());
            sliceMaxs = sliceDims.array() - 1;
        }
        
        /* Non-copyable */
        MatrixBaseXD(const MatrixBaseXD& mat) = delete;
        template<bool o, typename s, int rm, int cm>
        MatrixBaseXD(const MatrixBaseXD<o, s, rm, cm>& mat) = delete;
        MatrixBaseXD& operator=(const MatrixBaseXD& mat) = delete;
        
    private:
        mxArray*        mxData;
        Scalar*         matData;
        Dimensions      dimensions;
        Dimensions      sliceDims;
        Dimensions      sliceOffsets;
        Dimensions      sliceMins;
        Dimensions      sliceMaxs;
        int64_t      	numElements;
        int64_t       	sliceSize;
        
};

template<typename Scalar = double, int r = Eigen::Dynamic, int c = Eigen::Dynamic>
using MatrixXD = MatrixBaseXD<true, Scalar, r, c>;

template<typename Scalar = double, int r = Eigen::Dynamic, int c = Eigen::Dynamic>
using OutputMatrixXD = MatrixBaseXD<false, Scalar, r, c>;

/* ConstMatrix typedefs */
template<typename Scalar = double>
using ConstMatrixXDX = ConstMatrixXD<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/* Matrix typedefs */
template<typename Scalar = double>
using MatrixXDX = MatrixXD<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/* OutputMatrix typedefs */
template<typename Scalar = double>
using OutputMatrixXDX = OutputMatrixXD<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/* Double typedefs */
typedef ConstMatrixXDX<>        ConstMatrixXDXd;

typedef MatrixXDX<>             MatrixXDXd;

typedef OutputMatrixXDX<>    	OutputMatrixXDXd;

}   // namespace Mex

#endif
