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

#include "test_file_creation.h"

CPPUNIT_TEST_SUITE_REGISTRATION( TestFileCreation );

void TestFileCreation::setUp() {
}

void TestFileCreation::tearDown() {
}

void TestFileCreation::test_ply_file() {
    const auto ref = this->generate_mesh();
    this->assert_ply_shape("test.ply", ref);
}

void TestFileCreation::test_stl_file() {
    const auto ref = this->generate_mesh();
    this->assert_stl_shape("test.stl", ref);
}

void TestFileCreation::test_obj_file() {
    const auto ref = this->generate_mesh();
    this->assert_obj_shape("test.obj", ref);
}

TestFileCreation::MeshReference TestFileCreation::generate_mesh() const {
    // set number of threads to 1
    omp_set_num_threads(1);

    // read scalar field
    ScalarField sf("co_2pi_x.cub", ScalarFieldInputFileType::SFF_CUB);
    sf.read();

    // perform marching cubes algorithm
    IsoSurface is(&sf);
    is.marching_cubes(0.01);

    // construct mesh
    IsoSurfaceMesh ism(&sf, &is);
    ism.construct_mesh(true);

    // create files for all supported output formats
    ism.write_to_file("test.ply", "co molecule 2pi_x orbital test", "molecule");
    ism.write_to_file("test.stl", "co molecule 2pi_x orbital test", "molecule");
    ism.write_to_file("test.obj", "co molecule 2pi_x orbital test", "molecule");

    return {
        ism.get_vertices().size(),
        ism.get_normals().size(),
        ism.get_texcoords().size() / 6
    };
}

void TestFileCreation::assert_obj_shape(const std::string& filename, const MeshReference& ref) const {
    CPPUNIT_ASSERT(std::filesystem::exists(filename));

    std::ifstream infile(filename);
    CPPUNIT_ASSERT(infile.good());

    std::string line;
    size_t vertex_lines = 0;
    size_t normal_lines = 0;
    size_t face_lines = 0;

    while(std::getline(infile, line)) {
        if(line.rfind("v ", 0) == 0) {
            vertex_lines++;
        } else if(line.rfind("vn ", 0) == 0) {
            normal_lines++;
        } else if(line.rfind("f ", 0) == 0) {
            face_lines++;
        }
    }

    CPPUNIT_ASSERT_EQUAL(ref.vertices, vertex_lines);
    CPPUNIT_ASSERT_EQUAL(ref.normals, normal_lines);
    CPPUNIT_ASSERT_EQUAL(ref.triangles, face_lines);
}

void TestFileCreation::assert_ply_shape(const std::string& filename, const MeshReference& ref) const {
    CPPUNIT_ASSERT(std::filesystem::exists(filename));

    std::ifstream infile(filename, std::ios::binary);
    CPPUNIT_ASSERT(infile.good());

    std::string line;
    size_t header_size = 0;
    size_t header_vertices = 0;
    size_t header_faces = 0;
    bool seen_ply = false;
    bool seen_binary_format = false;

    while(std::getline(infile, line)) {
        header_size += line.size() + 1;

        if(line == "ply") {
            seen_ply = true;
        }
        if(line == "format binary_little_endian 1.0" || line == "format binary_big_endian 1.0") {
            seen_binary_format = true;
        }
        if(line.rfind("element vertex ", 0) == 0) {
            header_vertices = std::stoul(line.substr(15));
        }
        if(line.rfind("element face ", 0) == 0) {
            header_faces = std::stoul(line.substr(13));
        }
        if(line == "end_header") {
            break;
        }
    }

    CPPUNIT_ASSERT(seen_ply);
    CPPUNIT_ASSERT(seen_binary_format);
    CPPUNIT_ASSERT_EQUAL(ref.vertices, header_vertices);
    CPPUNIT_ASSERT_EQUAL(ref.triangles, header_faces);

    const size_t expected_binary_size =
        ref.vertices * (6 * sizeof(float)) +
        ref.triangles * (sizeof(uint8_t) + 3 * sizeof(unsigned int));

    CPPUNIT_ASSERT_EQUAL(header_size + expected_binary_size, (size_t)std::filesystem::file_size(filename));
}

void TestFileCreation::assert_stl_shape(const std::string& filename, const MeshReference& ref) const {
    CPPUNIT_ASSERT(std::filesystem::exists(filename));

    std::ifstream infile(filename, std::ios::binary);
    CPPUNIT_ASSERT(infile.good());

    // skip header
    infile.seekg(80, std::ios::beg);

    uint32_t triangle_count = 0;
    infile.read(reinterpret_cast<char*>(&triangle_count), sizeof(uint32_t));
    CPPUNIT_ASSERT(infile.good());

    CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(ref.triangles), triangle_count);

    const size_t expected_size = 80 + sizeof(uint32_t) + ref.triangles * (12 * sizeof(float) + sizeof(uint16_t));
    CPPUNIT_ASSERT_EQUAL(expected_size, (size_t)std::filesystem::file_size(filename));
}
