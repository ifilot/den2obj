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

#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>

#include "isosurface.h"

/**
 * @brief      structure to put glm::vec3 in a map
 */
struct KeyFuncs
{
    size_t operator()(const glm::vec3& k)const
    {
        return std::hash<float>()(k.x) ^ std::hash<float>()(k.y) ^ std::hash<float>()(k.z);
    }

    bool operator()(const glm::vec3& a, const glm::vec3& b)const
    {
            return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

/**
 * @brief      Class for iso surface mesh.
 */
class IsoSurfaceMesh{
private:
    std::unordered_map<glm::vec3, unsigned int, KeyFuncs, KeyFuncs> vertices_map;
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texcoords;
    std::vector<unsigned int> indices;

    const ScalarField* sf;
    const IsoSurface* is;

    glm::vec3 center;

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
     * @param[in]  center  whether to center structure
     */
    void construct_mesh(bool center);

    /**
     * @brief      write wavefront file
     *
     * @param[in]  filename  The filename
     * @param[in]  header    The header
     * @param[in]  name      The name
     */
    void write_obj(const std::string& filename, const std::string& header, const std::string& name);

private:
    /**
     * @brief      get the index of a vertex from unordered map
     *
     * @param[in]  v     vertex coordinates
     *
     * @return     the index
     */
    unsigned int get_index_vertex(const glm::vec3 v);

};

#endif //_OUTPUT_H
