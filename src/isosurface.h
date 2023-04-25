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

#ifndef _ISOSURFACE_H
#define _ISOSURFACE_H

#include <vector>
#include <iostream>
#include <cmath>
#include <mutex>

#include "edgetable.h"
#include "triangletable.h"
#include "scalar_field.h"
#include "math.h"

#define PRECISION_LIMIT 0.000000001

/*
 *            z
 *            |
 *            4----5
 *         /|     /|
 *        7-|--6 |
 *        | 0--|-1--- y
 *        |/     |/
 *        3----2
 *     /
 *    x
 *
 */

class Cube{
private:
    unsigned int i,j,k;                           // starting position of the cube
    float values[8];                              // vector holding cube values
    unsigned int cubidx;                          // cubeindex (whether there is an intersection)

public:
    Cube();

    Cube(unsigned int _i, unsigned int _j, unsigned int _k, const ScalarField &_vp);

    void set_cube_index(float _isovalue);

    unsigned int get_cube_index() const;

    float get_value_from_vertex(unsigned int _p) const;
    Vec3 get_position_from_vertex(unsigned int _p) const;
};

/*
 *               + 0
 *              /|\
 *             / | \
 *            /  |  \
 *           /   |   \
 *          /    |    \
 *         /     |     \
 *        +-------------+ 1
 *       3 \     |     /
 *          \    |    /
 *           \   |   /
 *            \  |  /
 *             \ | /
 *              \|/
 *               + 2
 */

class Tetrahedron{
private:
    unsigned int i,j,k;     //!< starting position of the cube
    float values[4];
    unsigned int tetidx;    //!< cubeindex (whether there is an intersection)
    Vec3 pos[4];       //!< coordinates of the gridpoints

public:
    Tetrahedron(unsigned int _i, unsigned int _j, unsigned int _k, const ScalarField &_vp, unsigned int _pos);

    void set_tetrahedron_index(float _isovalue);

    unsigned int get_tetrahedron_index() const;

    float get_value_from_vertex(unsigned int _p) const;

    const Vec3& get_position_from_vertex(unsigned int _p) const;
};

class Triangle{
public:
    Vec3 p1, p2, p3;

    Triangle();

    Triangle(const Vec3 &_p1, const Vec3 &_p2, const Vec3 &_p3);

    void transform_to_real(const ScalarField &_vp);

    float get_x(unsigned int i) const;

    float get_y(unsigned int i) const;

    float get_z(unsigned int i) const;
};

/**
 * @brief      generates an isosurface using either the marching cubes or the
 *             marching tetrahedra algorithm, input is a (tabulated) scalar
 *             field
 */
class IsoSurface {
private:
    std::vector<Cube> cube_table;
    std::vector<Tetrahedron> tetrahedra_table;
    std::vector<Triangle> triangles;
    ScalarField *vp_ptr;                // pointer to ScalarField obj
    std::array<unsigned int,3> grid_dimensions;
    float isovalue;                     // isovalue setting

public:
    /**
     * @brief      default constructor
     *
     * @param      _sf   pointer to ScalarField object
     */
    IsoSurface(ScalarField *_sf);

    /**
     * @brief      generate isosurface using marching cubes algorithm
     *
     * @param[in]  _isovalue  The isovalue
     */
    void marching_cubes(float _isovalue);

    /**
     * @brief      generate isosurface using marching tetrahedra algorithm
     *
     * @param[in]  _isovalue  The isovalue
     */
    void marching_tetrahedra(float _isovalue);

    /**
     * @brief      Gets the triangles pointer.
     *
     * @return     The triangles pointer.
     */
    const std::vector<Triangle>* get_triangles_ptr() const;

    inline float get_isovalue() const {
        return this->isovalue;
    }

    inline const ScalarField* get_scalar_field() const {
        return this->vp_ptr;
    }

private:
    void sample_grid_with_cubes(float _isovalue);
    void sample_grid_with_tetrahedra(float _isovalue);
    void construct_triangles_from_cubes(float _isovalue);
    void construct_triangles_from_tetrahedra(float _isovalue);
    Vec3 interpolate_from_cubes(const Cube &_cub, unsigned int _p1, unsigned int _p2, float _isovalue);
    Vec3 interpolate_from_tetrahedra(const Tetrahedron &_cub, unsigned int _p1, unsigned int _p2, float _isovalue);
};

#endif //_ISOSURFACE_H
