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

#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <fstream>
#include <set>
#include <vector>
#include <omp.h>

#include "math.h"
#include "isosurface.h"
#include "check_endian.h"

/**
 * @brief      structure to put Vec3 in a map
 */
struct KeyFuncs {
    size_t operator()(const Vec3& k) const {
        return std::hash<fpt>()(k[0]) ^ std::hash<fpt>()(k[1]) ^ std::hash<fpt>()(k[2]);
    }

    bool operator()(const Vec3& a, const Vec3& b) const {
            return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
    }
};

/**
 * @brief      Class for iso surface mesh.
 */
class IsoSurfaceMesh{
private:
    std::unordered_map<Vec3, unsigned int, KeyFuncs, KeyFuncs> vertices_map;
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<Vec2> texcoords;
    std::vector<unsigned int> indices;

    const ScalarField* sf;
    const IsoSurface* is;

    Vec3 center;

public:
    /**
     * @brief      build isosurface mesh object
     *
     * @param[in]  _sf   pointer to scalar field
     * @param[in]  _is   pointer to isosurface
     */
    IsoSurfaceMesh(const ScalarField* _sf, const IsoSurface* _is);

    /**
     * @brief      construct surface mesh
     *
     * @param[in]  center_mesh  whether to center structure
     */
    void construct_mesh(bool center_mesh);

    /**
     * @brief      write wavefront file
     *
     * @param[in]  filename  The filename
     * @param[in]  header    The header
     * @param[in]  name      The name
     */
    void write_obj(const std::string& filename, const std::string& header, const std::string& name);

    /**
     * @brief      write as binary ply file
     *
     * @param[in]  filename  The filename
     * @param[in]  header    The header
     * @param[in]  name      The name
     */
    void write_ply(const std::string& filename, const std::string& header, const std::string& name);

    /**
     * @brief      write as binary stl file
     *
     * @param[in]  filename  The filename
     */
    void write_stl(const std::string& filename);

    /**
     * @brief      get the vertices
     *
     * @return     vertices
     */
    inline const auto& get_vertices() const {
        return this->vertices;
    }

    /**
     * @brief      get the normals
     *
     * @return     normals
     */
    inline const auto& get_normals() const {
        return this->normals;
    }

    /**
     * @brief      get the texcoords
     *
     * @return     texcoords
     */
    inline const auto& get_texcoords() const {
        return this->texcoords;
    }

private:
    /**
     * @brief      get the index of a vertex from unordered map
     *
     * @param[in]  v     vertex coordinates
     *
     * @return     the index
     */
    unsigned int get_index_vertex(const Vec3 v);

};

#endif //_OUTPUT_H
