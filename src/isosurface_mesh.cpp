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

#include "isosurface_mesh.h"

/**
 * @brief      build isosurface mesh object
 *
 * @param[in]  _sf   pointer to scalar field
 * @param[in]  _is   pointer to isosurface
 */
IsoSurfaceMesh::IsoSurfaceMesh(const ScalarField* _sf,
                               const IsoSurface* _is) :
    sf(_sf),
    is(_is) {
}

/**
 * @brief      construct surface mesh
 *
 * @param[in]  center  whether to center structure
 */
void IsoSurfaceMesh::construct_mesh(bool center) {
   // grab center
    this->center = this->sf->get_mat_unitcell() * glm::vec3(0.5, 0.5, 0.5);

    for(unsigned int i=0; i<this->is->get_triangles_ptr()->size(); i++) {
        this->texcoords.push_back(glm::vec2(0,0));
        this->texcoords.push_back(glm::vec2(0,1));
        this->texcoords.push_back(glm::vec2(1,0));

        this->indices.push_back(this->get_index_vertex(is->get_triangles_ptr()->at(i).p1));
        this->indices.push_back(this->get_index_vertex(is->get_triangles_ptr()->at(i).p2));
        this->indices.push_back(this->get_index_vertex(is->get_triangles_ptr()->at(i).p3));
    }

    // build vertex vector from unordered map
    this->vertices.resize(this->vertices_map.size());
    for(auto it : this->vertices_map) {
        this->vertices[it.second] = it.first;
    }

    double dev = 0.01;
    this->normals.resize(this->vertices.size());

    // calculate normal vectors
    #pragma omp parallel for
    for(unsigned int i=0; i<this->vertices.size(); i++) {
        // get derivatives
        double dx0 = sf->get_value_interp(this->vertices[i][0] - dev, this->vertices[i][1], this->vertices[i][2]);
        double dx1 = sf->get_value_interp(this->vertices[i][0] + dev, this->vertices[i][1], this->vertices[i][2]);

        double dy0 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1] - dev, this->vertices[i][2]);
        double dy1 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1] + dev, this->vertices[i][2]);

        double dz0 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1], this->vertices[i][2] - dev);
        double dz1 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1], this->vertices[i][2] + dev);

        glm::vec3 normal((dx1 - dx0) / (2.0 * dev),
                    (dy1 - dy0) / (2.0 * dev),
                    (dz1 - dz0) / (2.0 * dev));
        normal = glm::normalize(normal);

        this->normals[i] = normal;
    }

    // center structure
    if(center) {
        float sx = 0.0f;
        float sy = 0.0f;
        float sz = 0.0f;

        #pragma omp parallel for reduction(+ : sx, sy, sz)
        for(unsigned int i=0; i<this->vertices.size(); i++) {
            sx += this->vertices[i][0];
            sy += this->vertices[i][1];
            sz += this->vertices[i][2];
        }

        glm::vec3 sum(sx, sy, sz);
        sum /= (float)this->vertices.size();

        #pragma omp parallel for
        for(unsigned int i=0; i<this->vertices.size(); i++) {
           this->vertices[i] -= sum;
        }
    }
}

/**
 * @brief      write wavefront file
 *
 * @param[in]  filename  The filename
 */
void IsoSurfaceMesh::write_obj(std::string filename) {
    std::cout << "Writing to " << filename << std::endl;
    std::ofstream myfile;
    myfile.open(filename.c_str());

    myfile << "# IsoTron OBJ File" << std::endl;
    myfile << "# www.ivofilot.nl" << std::endl;
    myfile << "o isosurface" << std::endl;

    for(unsigned int i=0; i<this->vertices.size(); i++) {
        myfile << "v " << this->vertices[i][0] << " " << this->vertices[i][1] << " " << this->vertices[i][2] << std::endl;
    }

    for(unsigned int i=0; i<this->normals.size(); i++) {
        myfile << "vn " << this->normals[i][0] << " " << this->normals[i][1] << " " << this->normals[i][2] << std::endl;
    }

    myfile << "s off" << std::endl;

    for(unsigned int i=0; i<this->indices.size(); i+=3) {
        myfile << "f " << this->indices[i]+1   << "//" << this->indices[i]+1
               << " "  << this->indices[i+1]+1 << "//" << this->indices[i+1]+1
               << " "  << this->indices[i+2]+1 << "//" << this->indices[i+2]+1 << std::endl;
    }

    myfile.close();
}

/**
 * @brief      get the index of a vertex from unordered map
 *
 * @param[in]  v     vertex coordinates
 *
 * @return     the index
 */
unsigned int IsoSurfaceMesh::get_index_vertex(const glm::vec3 v) {
    auto got = this->vertices_map.find(v);
    if(got != this->vertices_map.end()) {
        return got->second;
    } else {
        this->vertices_map.emplace(v, this->vertices_map.size());
        return this->get_index_vertex(v);
    }
}
