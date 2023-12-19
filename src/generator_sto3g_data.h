/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   DEN2OBJ is free software:                                            *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DEN2OBJ is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _GENERATOR_STO3G_DATA_H
#define _GENERATOR_STO3G_DATA_H

#include "math.h"

namespace GeneratorData {

static constexpr fpt basis_h[3][2] {
    {3.425251, 0.154329},
    {0.623914, 0.535328},
    {0.168855, 0.444635},
};

static constexpr fpt basis_c[3][3][2] {
    {
        {71.616837, 0.154329},
        {13.045096, 0.535328},
        {3.530512, 0.444635},
    },
    {
        {2.941249, -0.099967},
        {0.683483, 0.399513},
        {0.22229, 0.700115},
    },
    {
        {2.941249, 0.155916},
        {0.683483, 0.607684},
        {0.22229, 0.391957},
    }
};

} //  namespace GeneratorData

#endif // _GENERATOR_STO3G_DATA_H
