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
void IsoSurfaceMesh::construct_mesh(bool center_mesh) {
   // grab center
    this->center = this->sf->get_mat_unitcell() * Vec3(0.5, 0.5, 0.5);

    // build texture coordinates
    for(unsigned int i=0; i<this->is->get_triangles_ptr()->size(); i++) {
        this->texcoords.push_back(Vec2(0,0));
        this->texcoords.push_back(Vec2(0,1));
        this->texcoords.push_back(Vec2(1,0));

        // load all index vertices in a map; this operation needs to be done, else a SEGFAULT
        // will be thrown further down the lines
        this->get_index_vertex(is->get_triangles_ptr()->at(i).p1);
        this->get_index_vertex(is->get_triangles_ptr()->at(i).p2);
        this->get_index_vertex(is->get_triangles_ptr()->at(i).p3);
    }

    // check that 'get_index_vertex()' is being called for all vertices
    if(this->vertices_map.size() == 0 && this->is->get_triangles_ptr()->size() != 0) {
        throw std::runtime_error("Vertices map is empty, this is probably the result of the index vertices not being parsed. 'get_index_vertex() needs to be called for all vertices.");
    }

    // build vertex vector from unordered map
    this->vertices.resize(this->vertices_map.size());
    for(auto it : this->vertices_map) {
        this->vertices[it.second] = it.first;
    }

    double dev = 0.01;
    this->normals.resize(this->vertices.size());
    // std::cout << this->normals.size() << std::endl;

    // calculate normal vectors
    std::cout << "Calculating normal vectors using two-point stencil" << std::endl;
    #pragma omp parallel for schedule(static)
    for(unsigned int i=0; i<this->vertices.size(); i++) {
        // get derivatives
        double dx0 = sf->get_value_interp(this->vertices[i][0] - dev, this->vertices[i][1], this->vertices[i][2]);
        double dx1 = sf->get_value_interp(this->vertices[i][0] + dev, this->vertices[i][1], this->vertices[i][2]);

        double dy0 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1] - dev, this->vertices[i][2]);
        double dy1 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1] + dev, this->vertices[i][2]);

        double dz0 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1], this->vertices[i][2] - dev);
        double dz1 = sf->get_value_interp(this->vertices[i][0], this->vertices[i][1], this->vertices[i][2] + dev);

        Vec3 normal((dx1 - dx0) / (2.0 * dev),
                    (dy1 - dy0) / (2.0 * dev),
                    (dz1 - dz0) / (2.0 * dev));
        normal = -normal.normalized(); // the negative of the gradient is the correct normal

        this->normals[i] = normal * sgn(sf->get_value_interp(this->vertices[i][0], this->vertices[i][1], this->vertices[i][2]));
    }

    // build indices in right orientation based on face normal
    for(unsigned int i=0; i<this->is->get_triangles_ptr()->size(); i++) {
        // calculate face normal
        unsigned int id1 = this->get_index_vertex(is->get_triangles_ptr()->at(i).p1);
        unsigned int id2 = this->get_index_vertex(is->get_triangles_ptr()->at(i).p2);
        unsigned int id3 = this->get_index_vertex(is->get_triangles_ptr()->at(i).p3);

        // calculate the orientation of the face with respect to the normal
        const Vec3 face_normal = (this->normals[id1] + this->normals[id2] + this->normals[id3]) / 3.0f;
        const Vec3 orientation_face = ((this->vertices[id2] - this->vertices[id1]).cross(this->vertices[id3] - this->vertices[id1])).normalized();
        const float orientation = face_normal.dot(orientation_face);

        // if orientation is positive, the orientation is correct, if it is negative, the orientation is incorrect and two indices should be swapped
        if(orientation > 0.0f) {
            this->indices.push_back(id1);
            this->indices.push_back(id2);
            this->indices.push_back(id3);
        } else {
            this->indices.push_back(id2);
            this->indices.push_back(id1);
            this->indices.push_back(id3);
        }
    }

    // center structure if boolean is set
    if(center_mesh) {
        std::cout << "Centering structure" << std::endl;
        Vec3 sum = this->sf->get_mat_unitcell() * Vec3(0.5f, 0.5f, 0.5f);

        #pragma omp parallel for
        for(unsigned int i=0; i<this->vertices.size(); i++) {
           this->vertices[i] -= sum;
        }
    } else { // if not, use translation settings in Scalar Field (used in Cube files)

        // obtain translation values from ScalarField
        const Vec3 trans = this->sf->get_trans();

        #pragma omp parallel for
        for(unsigned int i=0; i<this->vertices.size(); i++) {
           this->vertices[i] += trans;
        }
    }
}

/**
 * @brief      Write to file, autorecognizes the file type
 *
 * @param[in]  filename  The filename
 * @param[in]  header    The header
 * @param[in]  name      The name
 */
void IsoSurfaceMesh::write_to_file(const std::string& filename, const std::string& header, const std::string& name) {
    boost::filesystem::path path(filename);

    if(filename.substr(filename.size()-4) == ".obj") {
        std::cout << "Writing mesh as Wavefront file (.obj)." << std::endl;
        this->write_obj(filename, path.filename().string(), path.filename().string());
    } else if(filename.substr(filename.size()-4) == ".ply") {
        std::cout << "Writing mesh as Standford Triangle Format file (.ply)." << std::endl;
        this->write_ply(filename, path.filename().string(), path.filename().string());
    } else if(filename.substr(filename.size()-4) == ".stl") {
        std::cout << "Writing mesh as Stereolithography file (.stl)." << std::endl;
        this->write_stl(filename);
    } else {
        throw std::runtime_error("Cannot interpret output file format. Please specify a valid extension.");
    }
}

/**
 * @brief      write wavefront file
 *
 * @param[in]  filename  The filename
 * @param[in]  header    The header
 * @param[in]  name      The name
 */
void IsoSurfaceMesh::write_obj(const std::string& filename, const std::string& header, const std::string& name) {
    std::ofstream outfile(filename);

    outfile << "# " << header << std::endl;
    outfile << "# Created by den2obj" << std::endl;
    outfile << "# https://github.com/ifilot/den2obj" << std::endl;
    outfile << "o " << name << std::endl;

    // calculate number of threads
    size_t nrthreads = omp_get_max_threads();
    omp_set_num_threads(nrthreads); // always allocate max threads
    std::stringstream local[nrthreads];

    // parallel writing vertices
    #pragma omp parallel
    {
        size_t threadnum = omp_get_thread_num();

        // calculate size
        size_t rem = this->vertices.size() % nrthreads;

        // divide task
        size_t start = this->vertices.size() / nrthreads * threadnum;
        size_t stop = this->vertices.size() / nrthreads * (threadnum + 1);
        if(threadnum == nrthreads - 1) {
            stop += rem;
        }

        char buffer[100];
        unsigned int cnt = 0;

        for(size_t i=start; i<stop; i++) {

            sprintf(buffer, "v %6.4f  %6.4f  %6.4f\n", this->vertices[i][0], this->vertices[i][1], this->vertices[i][2]);
            local[threadnum] << buffer;
        }
    }

    // merge results
    for(unsigned int i=0; i<nrthreads; i++) {
        outfile << local[i].str();
        local[i].str(std::string());    // clear stringstream
    }

    // parallel writing normals

    #pragma omp parallel
    {
        size_t threadnum = omp_get_thread_num();

        // calculate size
        size_t rem = this->normals.size() % nrthreads;

        // divide task
        size_t start = this->normals.size() / nrthreads * threadnum;
        size_t stop = this->normals.size() / nrthreads * (threadnum + 1);
        if(threadnum == nrthreads - 1) {
            stop += rem;
        }

        char buffer[100];
        unsigned int cnt = 0;

        for(size_t i=start; i<stop; i++) {

            sprintf(buffer, "vn %6.4f  %6.4f  %6.4f\n", this->normals[i][0], this->normals[i][1], this->normals[i][2]);
            local[threadnum] << buffer;
        }
    }

    // merge results
    for(unsigned int i=0; i<nrthreads; i++) {
        outfile << local[i].str();
        local[i].str(std::string());    // clear stringstream
    }

    outfile << "s off" << std::endl;

    // parallel writing faces

    #pragma omp parallel
    {
        size_t threadnum = omp_get_thread_num();

        // calculate size
        size_t rem = (this->indices.size() / 3) % nrthreads;

        // divide task
        size_t start = (this->indices.size() / 3) / nrthreads * threadnum * 3;
        size_t stop = (this->indices.size() / 3) / nrthreads * (threadnum+1) * 3;
        if(threadnum == nrthreads - 1) {
            stop += rem * 3;
        }

        char buffer[100];
        unsigned int cnt = 0;

        for(size_t i=start; i<stop; i+=3) {

            sprintf(buffer, "f %i//%i %i//%i %i//%i\n", this->indices[i]+1,   this->indices[i]+1,
                                                        this->indices[i+1]+1, this->indices[i+1]+1,
                                                        this->indices[i+2]+1, this->indices[i+2]+1);
            local[threadnum] << buffer;
        }
    }

    // merge results
    for(unsigned int i=0; i<nrthreads; i++) {
        outfile << local[i].str();
    }

    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    std::cout << "Writing as Wavefront (.obj) file: " << filename << " ("
    << (boost::format("%0.1f") % ((float)size / 1024.f)).str()
    << "kb)." << std::endl;
}

/**
 * @brief      write as binary ply file
 *
 * @param[in]  filename  The filename
 * @param[in]  header    The header
 * @param[in]  name      The name
 */
void IsoSurfaceMesh::write_ply(const std::string& filename, const std::string& header, const std::string& name) {
    std::ofstream outfile(filename, std::ios::binary);

    outfile << "ply" << std::endl;
    if(is_big_endian()) {
        outfile << "format binary_big_endian 1.0" << std::endl;
    } else {
        outfile << "format binary_little_endian 1.0" << std::endl;
    }

    outfile << "comment test" << std::endl;
    outfile << "element vertex " << this->vertices.size() << std::endl;
    outfile << "property float x" << std::endl;
    outfile << "property float y" << std::endl;
    outfile << "property float z" << std::endl;
    outfile << "property float nx" << std::endl;
    outfile << "property float ny" << std::endl;
    outfile << "property float nz" << std::endl;
    outfile << "element face " << (this->indices.size() / 3) << std::endl;
    outfile << "property list uchar uint vertex_indices" << std::endl;
    outfile << "end_header" << std::endl;

    // output vertex positions and normals
    for(unsigned int i=0; i<this->vertices.size(); i++) {
        outfile.write((char*)&this->vertices[i][0], sizeof(float) * 3);
        outfile.write((char*)&this->normals[i][0], sizeof(float) * 3);
    }

    // write indices
    static const uint8_t uchar_three = 3;
    for(unsigned int i=0; i<this->indices.size(); i+=3) {
        outfile.write((char*)&uchar_three, sizeof(uint8_t));
        outfile.write((char*)&this->indices[i], sizeof(unsigned int) * 3);
    }

    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    std::cout << "Writing as Stanford (.ply) file: " << filename << " ("
    << (boost::format("%0.1f") % ((float)size / 1024.f)).str()
    << "kb)." << std::endl;
}

/**
 * @brief      write as binary stl file
 *
 * @param[in]  filename  The filename
 */
void IsoSurfaceMesh::write_stl(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary);

    // for writing the header line
    static const char zerobyte = '\0';

    // write empty header of 80 bytes
    for(unsigned int i=0; i<80; i++) {
        outfile.write(&zerobyte, 1);
    }

    // write number of triangles
    uint32_t nr_triangles = this->indices.size() / 3;
    outfile.write((char*)&nr_triangles, sizeof(uint32_t));

    // for writing the attribute line
    static const uint16_t empty_attribute = 0;

    // loop over triangles
    for(unsigned int i=0; i<this->indices.size(); i+=3) {
        const unsigned int idx1 = this->indices[i];
        const unsigned int idx2 = this->indices[i+1];
        const unsigned int idx3 = this->indices[i+2];

        // The normals are based on the gradient of the underlying scalar
        // field and stored at the vertices. For STL files, the normal of
        // the faces is required. From the vertex normals, we take the sum
        // to get an approximation of the face normal. We calculate the normal
        // by taking the cross product of two triangle edges. Since the order
        // of the cross product determines the direction (sign), we multiply
        // that answer with the sign of the dot product with the sum of the
        // vertex normals. This will automatically flip the normal such that
        // it will point outwards.

        // calculate normal vector based on vertex normals
        Vec3 normal = (this->normals[idx1] +
                       this->normals[idx2] +
                       this->normals[idx3]);

        // calculate normal vector based on edges
        Vec3 edge1 = this->vertices[idx2] - this->vertices[idx1];
        Vec3 edge2 = this->vertices[idx3] - this->vertices[idx1];
        Vec3 normal_cross = (edge1.cross(edge2)).normalized();

        // flip outwards if necessary
        normal_cross *= sgn(normal.dot(normal_cross));

        // write normal vector
        outfile.write((char*)&normal_cross[0], 3 * sizeof(fpt));

        // write vertices
        outfile.write((char*)&vertices[idx1][0], 3 * sizeof(fpt));
        outfile.write((char*)&vertices[idx2][0], 3 * sizeof(fpt));
        outfile.write((char*)&vertices[idx3][0], 3 * sizeof(fpt));

        // write empty attribute value
        outfile.write((char*)&empty_attribute, sizeof(uint16_t));
    }

    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    std::cout << "Writing as Stereolithography (.stl) file: " << filename << " ("
    << (boost::format("%0.1f") % ((float)size / 1024.f)).str()
    << "kb)." << std::endl;
}

/**
 * @brief      get the index of a vertex from unordered map
 *
 * @param[in]  v     vertex coordinates
 *
 * @return     the index
 */
unsigned int IsoSurfaceMesh::get_index_vertex(const Vec3 v) {
    auto got = this->vertices_map.find(v);
    if(got != this->vertices_map.end()) {
        return got->second;
    } else {
        this->vertices_map.emplace(v, this->vertices_map.size());
        return this->get_index_vertex(v);
    }
}
