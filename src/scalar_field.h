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

#ifndef _SCALAR_FIELD_H
#define _SCALAR_FIELD_H

#include <string>
#include <vector>
#include <iostream>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <glm/glm.hpp>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include "float_parser.h"
#include "periodic_table.h"

template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

enum class OpenVDB_METHOD {
    ABSOLUTE,
    ABSOLUTE_LOG,
    POSITIVE,
    NEGATIVE,
    POSITIVE_LOG,
    NEGATIVE_LOG
};

class ScalarField{
private:
    std::string filename;
    float scalar;
    float mat[3][3];                //!< matrix dimensions
    float imat[3][3];               //!< inverse of matrix

    glm::mat3 mat33;                //!< glm version of the matrix
    glm::mat3 imat33;               //!< glm version of the inverse matrix

    float volume;                   //!< unit cell volume

    unsigned int grid_dimensions[3];
    std::vector<unsigned int> nrat;

    std::vector<glm::vec3> atom_pos;
    std::vector<unsigned int> atom_charges;

    std::vector<glm::vec3> atom_pos_exp;
    std::vector<unsigned int> atom_charges_exp;

    std::string gridline;
    std::vector<float> gridptr;         //!< grid to first pos of float array
    std::vector<float> gridptr2;        //!< grid to first pos of float array
    unsigned int gridsize;
    bool vasp5_input = false;
    bool has_read = false;
    bool header_read = false;
    std::ifstream infile;
    bool flag_is_locpot = false;         //!< whether scalar field is in LOCPOT style

    float trans[3];                      //!< translation vector for cube files

public:

    /**
     * @brief      constructor
     *
     * @param[in]  _filename        url to filename
     * @param[in]  _flag_is_locpot  whether this file is a locpot
     * @param[in]  _is_bin          is binary file
     */
    ScalarField(const std::string &_filename, bool _flag_is_locpot = false, bool _is_bin = false);

    /**
     * @brief      Gets the unitcell matrix.
     *
     * @return     The unitcell matrix.
     */
    glm::mat3 get_unitcell_matrix();

    /**
     * @brief      determines if this Scalarfield is in LOCPOT-style
     *
     * @return     True if locpot, False otherwise.
     */
    inline bool is_locpot() const {
        return this->flag_is_locpot;
    }

    void read();

    void read_header_and_atoms();

    /*
     * float get_value_interp(x,y,z)
     *
     * Grabs a value from the 3D scalar field. Calculate the value
     * by using a trilinear interpolation.
     *
     * The trilinear interpolation algorithm has been extracted from:
     * http://paulbourke.net/miscellaneous/interpolation/
     *
     * Future algorithm can make use of a cubic interpolation.
     *
     */
    float get_value_interp(float x, float y, float z) const;

    /**
     * @brief      test whether point is inside unit cell
     *
     * @param[in]  x     x position
     * @param[in]  y     y position
     * @param[in]  z     z position
     *
     * @return     True if inside, False otherwise.
     */
    bool is_inside(float x, float y, float z) const;

    float get_value(unsigned int i, unsigned int j, unsigned int k) const;

    glm::vec3 grid_to_realspace(float i, float j, float k) const;

    glm::vec3 realspace_to_grid(float i, float j, float k) const;

    glm::vec3 realspace_to_direct(float i, float j, float k) const;

    void copy_grid_dimensions(unsigned int _grid_dimensions[]) const;

    /**
     * @brief      Gets the maximum value in scalar field.
     *
     * @return     The maximum.
     */
    float get_max() const;

    /**
     * @brief      Gets the minimum value in scalar field.
     *
     * @return     The minimum.
     */
    float get_min() const;

    glm::vec3 get_atom_position(unsigned int atid) const;

    inline const glm::mat3& get_mat_unitcell() const {
        return this->mat33;
    }

    inline const glm::mat3& get_mat_unitcell_inverse() const {
        return this->imat33;
    }

    /**
     * @brief      Gets the grid pointer.
     *
     * @return     The grid pointer.
     */
    inline const float* get_grid_ptr() const {
        return &this->gridptr[0];
    }

    /**
     * @brief      Gets the number of data points in the grid
     *
     * @return     Number of data points
     */
    unsigned int get_size() const {
        return this->gridptr.size();
    }

    /**
     * @brief      Gets the filename.
     *
     * @return     The filename.
     */
    inline const std::string& get_filename() const {
        return this->filename;
    }

    /**
     * @brief      Gets the translation vector
     *
     * @return     Translation vector
     */
    inline const float* get_trans() const {
        return this->trans;
    }

    /**
     * @brief      Writes to an OpenVDB file
     *
     * @param[in]  filename  The filename
     * @param[in]  method    The method (absolute value, positive, negative, log)
     */
    void write_to_vdb(const std::string& filename, OpenVDB_METHOD method) const;

private:
    void test_vasp5();
    void read_scalar();
    void read_matrix();
    void read_grid_dimensions();
    void read_nr_atoms();
    void read_atom_positions();
    void read_grid();
    float get_max_direction(unsigned int dim);

    /**
     * @brief      Calculate the inverse of the unit cell matrix
     */
    void calculate_inverse();

    /**
     * @brief      Calculate volume of the unit cell
     */
    void calculate_volume();

    /**
     * @brief      Load a binary file
     */
    void load_binary();

    /**
     * @brief      Load a cube file
     */
    void load_cube_file();
};

#endif //_SCALAR_FIELD_H
