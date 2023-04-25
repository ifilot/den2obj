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

#include "test_isosurface.h"

CPPUNIT_TEST_SUITE_REGISTRATION( TestIsosurface );

void TestIsosurface::setUp() {
    this->sf = std::make_unique<ScalarField>("ch4.d2o", ScalarFieldInputFileType::SFF_D2O);
}

void TestIsosurface::tearDown() {
}

void TestIsosurface::test_marching_cubes() {
    // construct isosurface
    IsoSurface is(this->sf.get());
    is.marching_cubes(0.01);

    // construct mesh storage object
    IsoSurfaceMesh ism(sf.get(), &is);
    ism.construct_mesh(true);

    CPPUNIT_ASSERT_EQUAL((size_t)8170, ism.get_vertices().size());
    CPPUNIT_ASSERT_EQUAL((size_t)8170, ism.get_normals().size());
    CPPUNIT_ASSERT_EQUAL((size_t)48984, ism.get_texcoords().size());

    // get minimum and maximum value
    float minx = 0.0;
    float miny = 0.0;
    float minz = 0.0;
    float maxx = 0.0;
    float maxy = 0.0;
    float maxz = 0.0;
    for(unsigned int i=0; i<ism.get_vertices().size(); i++) {
        minx = std::min(minx, ism.get_vertices()[i][0]);
        miny = std::min(miny, ism.get_vertices()[i][1]);
        minz = std::min(minz, ism.get_vertices()[i][2]);

        maxx = std::max(maxx, ism.get_vertices()[i][0]);
        maxy = std::max(maxy, ism.get_vertices()[i][1]);
        maxz = std::max(maxz, ism.get_vertices()[i][2]);
    }

    // check minimum and maximum values; note that the ordering of the
    // vertices object is more-or-less random due to them being generated
    // under OpenMP
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.98314094543457, minx, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.98314094543457, miny, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.98314094543457, minz, 1e-8);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.98314094543457, maxx, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.98314094543457, maxy, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.98314094543457, maxz, 1e-8);

    // verify that all normals are, by approximation, of unit size
    for(unsigned int i=0; i<ism.get_vertices().size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, ism.get_normals()[i].norm(), 1e-6);
    }
}

void TestIsosurface::test_obj() {
    // construct isosurface
    IsoSurface is(this->sf.get());
    is.marching_cubes(0.01);

    // construct mesh storage object
    IsoSurfaceMesh ism(sf.get(), &is);
    ism.construct_mesh(true);

    std::string fname = "ch4.obj";
    if(std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    ism.write_obj(fname, "test", "test");
    CPPUNIT_ASSERT(std::filesystem::exists(fname));

    // .obj files can vary a little bit in size, so we have to sample for a range
    CPPUNIT_ASSERT_GREATEREQUAL((size_t)900000, std::filesystem::file_size(fname));
    CPPUNIT_ASSERT_LESSEQUAL((size_t)1000000, std::filesystem::file_size(fname));
}

void TestIsosurface::test_ply() {
    // construct isosurface
    IsoSurface is(this->sf.get());
    is.marching_cubes(0.01);

    // construct mesh storage object
    IsoSurfaceMesh ism(sf.get(), &is);
    ism.construct_mesh(true);

    std::string fname = "ch4.ply";
    if(std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    ism.write_ply(fname, "test", "test");
    CPPUNIT_ASSERT(std::filesystem::exists(fname));

    // because .ply files are written in binary form, their size is fixed
    CPPUNIT_ASSERT_EQUAL((size_t)408588, std::filesystem::file_size(fname));
}

void TestIsosurface::test_stl() {
    // construct isosurface
    IsoSurface is(this->sf.get());
    is.marching_cubes(0.01);

    // construct mesh storage object
    IsoSurfaceMesh ism(sf.get(), &is);
    ism.construct_mesh(true);

    std::string fname = "ch4.stl";
    if(std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    ism.write_stl(fname);
    CPPUNIT_ASSERT(std::filesystem::exists(fname));

    // because .stl files are written in binary form, their size is fixed
    CPPUNIT_ASSERT_EQUAL((size_t)816484, std::filesystem::file_size(fname));
}
