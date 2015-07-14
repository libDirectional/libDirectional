
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

#ifndef _MEX_TRAITS_H_
#define _MEX_TRAITS_H_

#include <mex.h>
#include <cstdint>

namespace Mex {

template<typename T>
struct Traits { };

/* Float */
template<>
struct Traits<float> {
    static const mxClassID MxClassID = mxSINGLE_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsSingle(array) && !mxIsSparse(array);
    }
};

/* Double */
template<>
struct Traits<double> {
    static const mxClassID MxClassID = mxDOUBLE_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsDouble(array) && !mxIsSparse(array);
    }
};

/* 8 bit integers */
template<>
struct Traits<int8_t> {
    static const mxClassID MxClassID = mxINT8_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsInt8(array) && !mxIsSparse(array);
    }
};

template<>
struct Traits<uint8_t> {
    static const mxClassID MxClassID = mxUINT8_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsUint8(array) && !mxIsSparse(array);
    }
};

/* 16 bit integers */
template<>
struct Traits<int16_t> {
    static const mxClassID MxClassID = mxINT16_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsInt16(array) && !mxIsSparse(array);
    }
};

template<>
struct Traits<uint16_t> {
    static const mxClassID MxClassID = mxUINT16_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsUint16(array) && !mxIsSparse(array);
    }
};

/* 32 bit integers */
template<>
struct Traits<int32_t> {
    static const mxClassID MxClassID = mxINT32_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsInt32(array) && !mxIsSparse(array);
    }
};

template<>
struct Traits<uint32_t> {
    static const mxClassID MxClassID = mxUINT32_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsUint32(array) && !mxIsSparse(array);
    }
};

/* 64 bit integers */
template<>
struct Traits<int64_t> {
    static const mxClassID MxClassID = mxINT64_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsInt64(array) && !mxIsSparse(array);
    }
};

template<>
struct Traits<uint64_t> {
    static const mxClassID MxClassID = mxUINT64_CLASS;
    
    static bool isValidArray(const mxArray* array) {
        return mxIsUint64(array) && !mxIsSparse(array);
    }
};

}   // namespace Mex

#endif
