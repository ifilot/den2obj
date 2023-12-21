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

    // create file
    ism.write_to_file("test.ply", "co molecule 2pi_x orbital test", "molecule");

    // test md5sum
    CPPUNIT_ASSERT_EQUAL(this->md5("test.ply"), std::string("ec143f6ebf9b72761803c3b091f54bf3"));
}

void TestFileCreation::test_stl_file() {
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

    // create file
    ism.write_to_file("test.stl", "co molecule 2pi_x orbital test", "molecule");

    // test md5sum
    CPPUNIT_ASSERT_EQUAL(this->md5("test.stl"), std::string("c2194ba639caf5092654862bb9f93298"));
}

void TestFileCreation::test_obj_file() {
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

    // create file
    ism.write_to_file("test.obj", "co molecule 2pi_x orbital test", "molecule");

    // test md5sum
    CPPUNIT_ASSERT_EQUAL(this->md5("test.obj"), std::string("e2b3e09f9c010dac99a7bc0137c187ec"));
}

/**
 * @brief      Calculate MD5 checksum of a file
 *
 * @param[in]  name      Path to the file
 * 
 * @return     32 byte string containing md5 checksum
 */
std::string TestFileCreation::md5(const std::string& filename) {
    // read the file
    std::ifstream mfile(filename, std::ios::binary | std::ios::ate);
    std::streamsize size = mfile.tellg();
    mfile.seekg(0, std::ios::beg);
    char buffer[size];
    mfile.read(buffer, size);
    mfile.close();

    // output variable for hash
    unsigned char hash[MD5_DIGEST_LENGTH];

    // calculate the md5 hash
    EVP_MD_CTX *ctx = EVP_MD_CTX_create();
    EVP_MD_CTX_init(ctx);
    const EVP_MD *md_type = EVP_md5();
    EVP_DigestInit_ex(ctx, md_type, NULL);
    EVP_DigestUpdate(ctx, buffer, size);
    EVP_DigestFinal_ex(ctx, hash, NULL);
    EVP_MD_CTX_destroy(ctx);

    // output as a 32-byte hex-string
    std::stringstream ss;
    for(int i = 0; i < MD5_DIGEST_LENGTH; i++){
        ss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>( hash[i] );
    }
    return ss.str();
}