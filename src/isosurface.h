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

#define PRECISION_LIMIT 0.000000001

/* each cube holds 8 pointers to values (at the vertices) and the starting
location defined as the lowest i, j and k value for the cube */

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
    /* cube constructor takes the grid coordinate for the cube origin and extracts
    the values for the vertices from the grids file */
    Cube(unsigned int _i, unsigned int _j, unsigned int _k,
             const ScalarField &_vp);
    /* set the cubidx value by comparing the value at the vertices with the
    supplied isovalue */
    void set_cube_index(float _isovalue);

    unsigned int get_cube_index() const;

    float get_value_from_vertex(unsigned int _p) const;
    glm::vec3 get_position_from_vertex(unsigned int _p) const;
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
    unsigned int i,j,k;                                // starting position of the cube
    float values[4];
    unsigned int tetidx;                             // cubeindex (whether there is an intersection)
    glm::vec3 pos[4];                                        // coordinates of the gridpoints

public:
    /* tetrahedron constructor takes the grid coordinate for the tetrahedron
    origin and extracts the values for the vertices from the grids file */
    Tetrahedron(unsigned int _i, unsigned int _j, unsigned int _k,
             const ScalarField &_vp, unsigned int _pos);
    /* set the cubidx value by comparing the value at the vertices with the
    supplied isovalue */
    void set_tetrahedron_index(float _isovalue);

    unsigned int get_tetrahedron_index() const;

    float get_value_from_vertex(unsigned int _p) const;
    const glm::vec3& get_position_from_vertex(unsigned int _p) const;
};

class Triangle{
public:
    glm::vec3 p1, p2, p3;

    Triangle();
    Triangle(const glm::vec3 &_p1, const glm::vec3 &_p2, const glm::vec3 &_p3);
    void transform_to_real(const ScalarField &_vp);
    float get_x(unsigned int i) const;
    float get_y(unsigned int i) const;
    float get_z(unsigned int i) const;
};

/* generates an isosurface using either the marching cubes or the marching
tetrahedra algorithm, input is a (tabulated) scalar field */

class IsoSurface {
private:
    std::vector<Cube> cube_table;
    std::vector<Tetrahedron> tetrahedra_table;
    std::vector<Triangle> triangles;
    ScalarField *vp_ptr;                // pointer to ScalarField obj
    unsigned int grid_dimensions[3];
    float isovalue;                     // isovalue setting

public:
    IsoSurface(ScalarField *_sf);                           // default constructor
    void marching_cubes(float _isovalue);                   // generate isosurface
    void marching_tetrahedra(float _isovalue);              // generate isosurface
    const std::vector<Triangle>* get_triangles_ptr() const; // return pointer to the triangles vector

    inline float get_isovalue() const {
        return this->isovalue;
    }

    inline const ScalarField* get_scalar_field() const {
        return this->vp_ptr;
    }

private:
    /* sample the grid and collect all voxels (cubes) where the cubidx
    value is not 0 or 256 (i.e. all voxels which are intersected by the isovalue */
    void sample_grid_with_cubes(float _isovalue);
    void sample_grid_with_tetrahedra(float _isovalue);
    /* extract the triangles from the list of voxels */
    void construct_triangles_from_cubes(float _isovalue);
    void construct_triangles_from_tetrahedra(float _isovalue);
    /* construct the glm::vec3 coordinate of the triangle */
    glm::vec3 interpolate_from_cubes(const Cube &_cub, unsigned int _p1,
        unsigned int _p2, float _isovalue);
    glm::vec3 interpolate_from_tetrahedra(const Tetrahedron &_cub, unsigned int _p1,
        unsigned int _p2, float _isovalue);
};

#endif //_ISOSURFACE_H
