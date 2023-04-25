/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   EDP is free software:                                                *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   EDP is distributed in the hope that it will be useful,               *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _MATRICES_H
#define _MATRICES_H

#include <Eigen/Dense>

typedef float fpt;  // general floating point type
typedef Eigen::Matrix<fpt, 3, 3, Eigen::RowMajor> MatrixUnitcell;

typedef Eigen::Matrix<fpt, 3, 1> Vec3;
typedef Eigen::Matrix<fpt, 2, 1> Vec2;

template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif // _MATRICES_H
