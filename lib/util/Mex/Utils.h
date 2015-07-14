
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

#ifndef _MEX_UTILS_H_
#define _MEX_UTILS_H_

#include "../../external/Eigen/Dense"
#include <stdexcept>

namespace Mex {

typedef Eigen::Matrix<int64_t, Eigen::Dynamic, 1> Dimensions;

struct Utils {
    static Dimensions getDimensions(const mxArray* array) {
        mxAssert(!mxIsSparse(array), "Array must be dense.");
        
        const mwSize numDims = mxGetNumberOfDimensions(array);
        const mwSize* dims   = mxGetDimensions(array);
        
        Dimensions outDims(numDims);
        
        for (mwSize i = 0; i < numDims; ++i) {
            outDims(i) = (int64_t)dims[i];
        }
        
        return outDims;
    }
    
    // from the slice size vectors [..., i, ...] and [..., j, ...]
    // creates the slice size vector [..., max(i, j), ...], taking 
    // into account that the inputs might have different sizes
    static Dimensions expandSliceDims(const Dimensions& sliceDimsA,
                                      const Dimensions& sliceDimsB) {
        Dimensions sliceDims(std::max(sliceDimsA.size(), sliceDimsB.size()));
        
        for (int64_t i = 0; i < sliceDims.size(); i++) {
            sliceDims(i) = std::max(i < sliceDimsA.size() ? sliceDimsA(i) : 1, 
                                    i < sliceDimsB.size() ? sliceDimsB(i) : 1);
        }
        
        return sliceDims;
    }
    
    template<class MatA, class MatB>
    static Dimensions expandSlices(const MatA& matA,
                                   const MatB& matB) {
        return expandSliceDims(matA.slices(), matB.slices());
    }
    
    static Dimensions expandSlices(const Dimensions& dimsA,
                                   const Dimensions& dimsB) {
        return expandSliceDims(dimsA.tail(dimsA.size() - 2),
                               dimsB.tail(dimsB.size() - 2));
    }
    
    // checks if sliceMins <= slice <= sliceMax coefficient wise,
    // taking into account that the inputs might have different sizes
    static bool isValidSlice(const Dimensions& slice,
                             const Dimensions& sliceMins,
                             const Dimensions& sliceMaxs) {
        const Dimensions::Index sliceDims = sliceMins.size();
        const Dimensions sHead = slice.head(sliceDims);
        const Dimensions sTail = slice.tail(slice.size() - sliceDims);
        
        if (slice.size() < sliceDims) {
            return false;
        } else if (sHead.size() == 0) {
            return sTail.isZero();
        } else {
            return sliceMins.cwiseMin(sHead) == sliceMins &&
                   sliceMaxs.cwiseMax(sHead) == sliceMaxs &&
                   (slice.size() == sliceDims || sTail.isZero());
        }
    }
    
    template<int r>
    static void checkRows(int64_t rows) {
        if (r != Eigen::Dynamic && r != rows) {
            throw std::invalid_argument("Mismatch between given (" +
                                        std::to_string(rows) +
                                        ") and expected (" +
                                        std::to_string(r) +
                                        ") number of rows.");
        }
    }
    
    template<int c>
    static void checkCols(int64_t cols) {
        if (c != Eigen::Dynamic && c != cols) {
            throw std::invalid_argument("Mismatch between given (" +
                                        std::to_string(cols) +
                                        ") and expected (" +
                                        std::to_string(c) +
                                        ") number of columns.");
        }
    }
    
    template<typename Scalar, int r, int c, bool isMatrix = true>
    static Scalar* checkArray(mxArray* array) {
        if (!Traits<Scalar>::isValidArray(array)) {
            throw std::invalid_argument("MX array of invalid type.");
        }
        
        if (isMatrix && mxGetNumberOfDimensions(array) != 2) {
            throw std::invalid_argument("MX array is not a matrix.");
        }
        
        const mwSize* dims = mxGetDimensions(array);
        
        checkRows<r>((int64_t)dims[0]);
        checkCols<c>((int64_t)dims[1]);
        
        return (Scalar*) mxGetData(array);
    }
    
    template<typename Scalar, int r, int c, bool isMatrix = true>
    static const Scalar* checkArray(const mxArray* array) {
        if (!Traits<Scalar>::isValidArray(array)) {
            throw std::invalid_argument("MX array of invalid type.");
        }
        
        if (isMatrix && mxGetNumberOfDimensions(array) != 2) {
            throw std::invalid_argument("MX array is not a matrix.");
        }
        
        const mwSize* dims = mxGetDimensions(array);
        
        checkRows<r>((int64_t)dims[0]);
        checkCols<c>((int64_t)dims[1]);
        
        return (const Scalar*) mxGetData(array);
    }

};

}   // namespace Mex

#endif
