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

#include "scalar_field.h"

/**
 * @brief      constructor
 *
 * @param[in]  _filename   url to filename
 * @param[in]  _flag_is_locpot  whether this file is a locpot
 */
ScalarField::ScalarField(const std::string &_filename, bool _flag_is_locpot, bool _is_bin) {
    this->filename = _filename;
    this->trans[0] = 0.0;
    this->trans[1] = 0.0;
    this->trans[2] = 0.0;

    if (!boost::filesystem::exists(this->filename)) {
        throw std::runtime_error("Cannot open " + this->filename + "!");
    }

    if(this->filename.substr(this->filename.size()-4) == ".cub") {
        this->load_cube_file();
    } else if(_is_bin) {
        this->load_binary();
    } else {
        this->scalar = -1;
        this->vasp5_input = false;
        this->has_read = false;
        this->header_read = false;
        this->flag_is_locpot = _flag_is_locpot;
    }
}

/**
 * @brief      Gets the unitcell matrix.
 *
 * @return     The unitcell matrix.
 */
glm::mat3 ScalarField::get_unitcell_matrix() {
    glm::mat3 out;
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            out[i][j] = mat[i][j];
        }
    }

    return out;
}

/*
 * void read()
 *
 * Wrapper function that reads in the OUTCAR file
 *
 * Usage: sf.read(true);
 *
 */
void ScalarField::read_header_and_atoms() {
    if(this->header_read) {
        return;
    }

    this->test_vasp5();
    this->read_scalar();
    this->read_matrix();
    this->read_nr_atoms();
    this->read_atom_positions();
    this->read_grid_dimensions();
}

/*
 * void read()
 *
 * Wrapper function that reads in the OUTCAR file
 *
 * Usage: sf.read(true);
 *
 */
void ScalarField::read() {
    if(this->has_read) {
        return;
    }

    this->read_header_and_atoms();
    this->read_grid();
}

/*
 * void test_vasp5()
 *
 * Test if the input file is a VASP5 output file
 *
 */
void ScalarField::test_vasp5() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    for(unsigned int i=0; i<5; i++) { // discard first two lines
        std::getline(infile, line);
    }
    std::getline(infile, line);
    // check if this line contains atomic information (i.e. alpha-characters)
    boost::regex regex_vasp_version("^(.*[A-Za-z]+.*)$");
    boost::smatch what;
    if(boost::regex_match(line, what, regex_vasp_version)) {
        this->vasp5_input = true;
    }
}

/*
 * void read_scalar()
 *
 * Read the scalar value from the 2nd line of the
 * CHGCAR file. Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_scalar() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    std::getline(infile, line); // discard this line

    std::getline(infile, line);
    boost::regex regex_scalar("^\\s*([0-9.-]+)\\s*$");
    boost::smatch what;
    if(boost::regex_match(line, what, regex_scalar)) {
        this->scalar = boost::lexical_cast<float>(what[1]);
    } else {
        this->scalar = -1;
    }
}

/*
 * void read_matrix()
 *
 * Reads the matrix that defines the unit cell
 * in the CHGCAR file. The inverse of that matrix
 * is automatically constructed.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_matrix() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    for(unsigned int i=0; i<2; i++) { // discard first two lines
        std::getline(infile, line);
    }

    // setup match pattern
    boost::regex regex_vasp_matrix_line("^\\s*([0-9.-]+)\\s+([0-9.-]+)\\s+([0-9.-]+)\\s*$");
    for(unsigned int i=0; i<3; i++) {
        std::getline(infile, line);
        boost::smatch what;
        if(boost::regex_match(line, what, regex_vasp_matrix_line)) {
            for(unsigned int j=0; j<3; j++) {
                mat[i][j] = boost::lexical_cast<float>(what[j+1]);
            }
        }
    }

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            this->mat[i][j] *= this->scalar;
        }
    }

    // also construct inverse matrix
    this->calculate_inverse();

    // calculate matrix volume
    this->calculate_volume();
}

/*
 * void read_nr_atoms()
 *
 * Read the number of atoms of each element. These
 * numbers are used to skip the required amount of
 * lines.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_nr_atoms() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    for(unsigned int i=0; i<5; i++) {
        std::getline(infile, line);
    }

    // store atom types
    if(this->vasp5_input) {
        std::getline(infile, line);
        std::vector<std::string> pieces;
        boost::trim(line);
        boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(const auto piece: pieces) {
            this->atom_charges.push_back(PeriodicTable::get().get_elnr(piece));
        }
    }

    // skip another line
    std::getline(infile, line);

    std::vector<std::string> pieces;
    boost::trim(line);
    boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(unsigned int i=0; i<pieces.size(); i++) {
        boost::trim(pieces[i]);
        this->nrat.push_back(boost::lexical_cast<unsigned int>(pieces[i]));
    }

    infile.close();
}

void ScalarField::read_atom_positions() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    // skip lines that contain atoms
    for(unsigned int i=0; i<(this->vasp5_input ? 8 : 7); i++) {
        std::getline(infile, line);
    }

    // read the atom positions
    for(unsigned int i=0; i<this->nrat.size(); i++) {
        for(unsigned int j=0; j<this->nrat[i]; j++) {
                std::getline(infile, line);
                std::vector<std::string> pieces;
                boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
                this->atom_pos.push_back(glm::vec3(boost::lexical_cast<float>(pieces[1]), boost::lexical_cast<float>(pieces[2]), boost::lexical_cast<float>(pieces[3])));
        }
    }

    infile.close();
}

/*
 * void read_grid_dimensions()
 *
 * Read the number of gridpoints in each
 * direction.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_grid_dimensions() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    // skip lines
    for(unsigned int i=0; i<(this->vasp5_input ? 10 : 9); i++) {
        std::getline(infile, line);
    }

    // // skip atom positions
    for(unsigned int i=0; i<this->nrat.size(); i++) {
        for(unsigned int j=0; j<this->nrat[i]; j++) {
            std::getline(infile, line);
        }
    }

    boost::trim(line);
    this->gridline = line;

    std::vector<std::string> pieces;
    boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(unsigned int i=0; i<pieces.size(); i++) {
        this->grid_dimensions[i] = boost::lexical_cast<unsigned int>(pieces[i]);
    }

    this->gridsize = this->grid_dimensions[0] * this->grid_dimensions[1] * this->grid_dimensions[2];
    this->header_read = true;

    infile.close();
}

/*
 * void read_grid()
 *
 * Read all the grid points. This function depends
 * on the the gridsize being set via the
 * read_grid_dimensions() function.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_grid() {
    this->read_header_and_atoms();

    this->infile.open(this->filename.c_str());
    std::string line;
    // skip irrelevant lines
    for(unsigned int i=0; i<(this->vasp5_input ? 10 : 9); i++) {
        std::getline(this->infile, line);
    }
    for(unsigned int i=0; i<this->atom_pos.size(); i++) {
        std::getline(this->infile, line);
    }

    float_parser p;

    /* read spin up */
    unsigned int linecounter=0; // for the counter
    static const boost::regex regex_augmentation("augmentation.*");

    while(std::getline(this->infile, line)) {
        // stop looping when a second gridline appears (this
        // is where the spin down part starts)
        if(line.compare(this->gridline) == 0) {
            std::cout << "I am breaking the loop" << std::endl;
            break;
        }

        boost::smatch what;
        if(boost::regex_match(line, what, regex_augmentation)) {
            std::cout << "Augmentation break encountered" << std::endl;
            break;
        }

        // set iterators
        std::string::const_iterator b = line.begin();
        std::string::const_iterator e = line.end();

        // parse
        std::vector<float> floats;
        boost::spirit::qi::phrase_parse(b, e, p, boost::spirit::ascii::space, floats);

        // expand gridptr with the new size
        unsigned int cursize = this->gridptr.size();
        this->gridptr.resize(cursize + floats.size());

        // For CHGCAR type files, the electron density is multiplied by the cell volume
        // as described by the link below:
        // https://cms.mpi.univie.ac.at/vasp/vasp/CHGCAR_file.html#file-chgcar
        // Hence, for these files, we have to divide the value at the grid point by the
        // cell volume. For LOCPOT files, we should *not* do this.
        if(this->flag_is_locpot) {      // LOCPOT type files
            for(unsigned int j=0; j<floats.size(); j++) {
                this->gridptr[cursize + j] = floats[j];
            }
        } else {    // CHGCAR type files
            for(unsigned int j=0; j<floats.size(); j++) {
                this->gridptr[cursize + j] = floats[j] / this->volume;
            }
        }

        linecounter++;

        if(this->gridptr.size() >= this->gridsize) {
            this->has_read = true;
        }
    }

    infile.close();
}

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
float ScalarField::get_value_interp(float x, float y, float z) const {
    if(!this->is_inside(x,y,z)) {
        return 0.0f;
    }

    // cast the input to the nearest grid point
    glm::vec3 r = this->realspace_to_grid(x,y,z);
    glm::vec3 d = this->realspace_to_direct(x,y,z);

    // calculate value using trilinear interpolation
    float xd = remainderf(r[0], 1.0);
    float yd = remainderf(r[1], 1.0);
    float zd = remainderf(r[2], 1.0);

    if(xd < 0) xd += 1.0;
    if(yd < 0) yd += 1.0;
    if(zd < 0) zd += 1.0;

    float x0 = floor(r[0]);
    float x1 = ceil(r[0]);
    float y0 = floor(r[1]);
    float y1 = ceil(r[1]);
    float z0 = floor(r[2]);
    float z1 = ceil(r[2]);

    return
    this->get_value(x0, y0, z0) * (1.0 - xd) * (1.0 - yd) * (1.0 - zd) +
    this->get_value(x1, y0, z0) * xd                 * (1.0 - yd) * (1.0 - zd) +
    this->get_value(x0, y1, z0) * (1.0 - xd) * yd                 * (1.0 - zd) +
    this->get_value(x0, y0, z1) * (1.0 - xd) * (1.0 - yd) * zd                 +
    this->get_value(x1, y0, z1) * xd                 * (1.0 - yd) * zd                 +
    this->get_value(x0, y1, z1) * (1.0 - xd) * yd                 * zd                 +
    this->get_value(x1, y1, z0) * xd                 * yd                 * (1.0 - zd) +
    this->get_value(x1, y1, z1) * xd                 * yd                 * zd;
}

/**
 * @brief      test whether point is inside unit cell
 *
 * @param[in]  x     x position
 * @param[in]  y     y position
 * @param[in]  z     z position
 *
 * @return     True if inside, False otherwise.
 */
bool ScalarField::is_inside(float x, float y, float z) const {
    // cast the input to the nearest grid point
    glm::vec3 d = this->realspace_to_direct(x,y,z);

    if(d[0] < 0 || d[0] > 1.0) {
        return false;
    }
    if(d[1] < 0 || d[1] > 1.0) {
        return false;
    }
    if(d[2] < 0 || d[2] > 1.0) {
        return false;
    }

    return true;
}

/*
 * float get_max_direction(dim)
 *
 * Get the maximum value in a particular dimension. This is a convenience
 * function for the get_value_interp() function.
 *
 */
float ScalarField::get_max_direction(unsigned int dim) {
    float sum = 0;
    for(unsigned int i=0; i<3; i++) {
        sum += this->mat[i][dim];
    }
    return sum;
}

/**
 * @brief      Calculate the inverse of the unit cell matrix
 */
void ScalarField::calculate_inverse() {
    float det = 0;
    for(unsigned int i=0;i<3;i++) {
        det += (this->mat[0][i]*(this->mat[1][(i+1)%3]*this->mat[2][(i+2)%3] - this->mat[1][(i+2)%3]*this->mat[2][(i+1)%3]));
    }

    for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++) {
                     this->imat[i][j] = ((this->mat[(i+1)%3][(j+1)%3] * this->mat[(i+2)%3][(j+2)%3]) - (this->mat[(i+1)%3][(j+2)%3]*this->mat[(i+2)%3][(j+1)%3]))/ det;
            }
     }
}

/**
 * @brief      Calculate volume of the unit cell
 */
void ScalarField::calculate_volume() {
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            this->mat33[i][j] = this->mat[i][j];
        }
    }

    this->imat33 = glm::inverse(mat33);

    this->volume = glm::dot(glm::cross(this->mat33[0], this->mat33[1]), this->mat33[2]);
}

/*
 * float get_value(i,j,k)
 *
 * Grabs the value at a particular grid point.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
float ScalarField::get_value(unsigned int i, unsigned int j, unsigned int k) const {
    unsigned int idx = k * this->grid_dimensions[0] * this->grid_dimensions[1] +
                       j * this->grid_dimensions[0] +
                       i;
    return this->gridptr[idx];
}

/*
 * glm::vec3 grid_to_realspace(i,j,k)
 *
 * Converts a grid point to a realspace vector. This function
 * is not being used at the moment.
 *
 */
glm::vec3 ScalarField::grid_to_realspace(float i, float j, float k) const {
    float dx = (float)i / (float)grid_dimensions[0];
    float dy = (float)j / (float)grid_dimensions[1];
    float dz = (float)k / (float)grid_dimensions[2];

    glm::vec3 r;
    r[0] = mat[0][0] * dx + mat[1][0] * dy + mat[2][0] * dz;
    r[1] = mat[0][1] * dx + mat[1][1] * dy + mat[2][1] * dz;
    r[2] = mat[0][2] * dx + mat[1][2] * dy + mat[2][2] * dz;

    return r;
}

/*
 * glm::vec3 realspace_to_grid(i,j,k)
 *
 * Convert 3d realspace vector to a position on the grid. Non-integer
 * values (i.e. floating point) are given as the result.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
glm::vec3 ScalarField::realspace_to_grid(float i, float j, float k) const {
    glm::vec3 r;
    r[0] = imat[0][0] * i + imat[0][1] * j + imat[0][2] * k;
    r[1] = imat[1][0] * i + imat[1][1] * j + imat[1][2] * k;
    r[2] = imat[2][0] * i + imat[2][1] * j + imat[2][2] * k;

    r[0] *= float(this->grid_dimensions[0]-1);
    r[1] *= float(this->grid_dimensions[1]-1);
    r[2] *= float(this->grid_dimensions[2]-1);

    return r;
}

/*
 * glm::vec3 realspace_to_direct(i,j,k)
 *
 * Convert 3d realspace vector to direct position.
 *
 */
glm::vec3 ScalarField::realspace_to_direct(float i, float j, float k) const {
    glm::vec3 r;
    r[0] = imat[0][0] * i + imat[0][1] * j + imat[0][2] * k;
    r[1] = imat[1][0] * i + imat[1][1] * j + imat[1][2] * k;
    r[2] = imat[2][0] * i + imat[2][1] * j + imat[2][2] * k;

    return r;
}

void ScalarField::copy_grid_dimensions(unsigned int _grid_dimensions[]) const {
    for(unsigned int i=0; i<3; i++) {
        _grid_dimensions[i] = this->grid_dimensions[i];
    }
}

/**
 * @brief      Gets the maximum value in scalar field.
 *
 * @return     The maximum.
 */
float ScalarField::get_max() const {
    return *std::max_element(this->gridptr.begin(), this->gridptr.end());
}

/**
 * @brief      Gets the minimum value in scalar field.
 *
 * @return     The minimum.
 */
float ScalarField::get_min() const {
    return *std::min_element(this->gridptr.begin(), this->gridptr.end());
}

glm::vec3 ScalarField::get_atom_position(unsigned int atid) const {
    if(atid < this->atom_pos.size()) {
        return this->mat33 * this->atom_pos[atid];
    } else {
        throw std::runtime_error("Requested atom id lies outside bounds");
    }
}

/**
 * @brief      Writes to an OpenVDB file
 *
 * @param[in]  filename  The filename
 * @param[in]  method    The method (absolute value, positive, negative, log)
 */
void ScalarField::write_to_vdb(const std::string& filename, OpenVDB_METHOD method) const {
    openvdb::initialize();

    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();

    const openvdb::Vec3f c(0.0f, 0.0f, 0.0f);

    openvdb::Coord ijk;
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

    for(unsigned int i=0; i<this->grid_dimensions[0]; i++) {    // x
        ijk[0] = i - this->grid_dimensions[0] / 2;
        for(unsigned int j=0; j<this->grid_dimensions[1]; j++) {    // y
            ijk[1] = j - this->grid_dimensions[1] / 2;
            for(unsigned int k=0; k<this->grid_dimensions[2]; k++) {    // z
                ijk[2] = k - this->grid_dimensions[2] / 2;
                unsigned int idx = k * this->grid_dimensions[0] * this->grid_dimensions[1] + j * this->grid_dimensions[0] + i;

                switch(method) {
                    case OpenVDB_METHOD::ABSOLUTE:
                        accessor.setValue(ijk, this->gridptr[idx] * this->gridptr[idx]);
                    break;
                    case OpenVDB_METHOD::ABSOLUTE_LOG:
                        accessor.setValue(ijk, std::log(this->gridptr[idx] * this->gridptr[idx]));
                    break;
                    case OpenVDB_METHOD::POSITIVE:
                        accessor.setValue(ijk, std::max(0.0f, sgn(this->gridptr[idx]) * this->gridptr[idx] * this->gridptr[idx]));
                    break;
                    case OpenVDB_METHOD::NEGATIVE:
                        accessor.setValue(ijk, std::fabs(std::min(0.0f, sgn(this->gridptr[idx]) * this->gridptr[idx] * this->gridptr[idx])));
                    break;
                    case OpenVDB_METHOD::POSITIVE_LOG:
                        accessor.setValue(ijk, std::log(std::max(0.0f, sgn(this->gridptr[idx]) * this->gridptr[idx] * this->gridptr[idx])));
                    break;
                    case OpenVDB_METHOD::NEGATIVE_LOG:
                        accessor.setValue(ijk, std::log(std::fabs(std::min(0.0f, sgn(this->gridptr[idx]) * this->gridptr[idx] * this->gridptr[idx]))));
                    break;
                }
            }
        }
    }

    // Identify the grid as a level set.
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    grid->setName("density");

    // Create a VDB file object and write out the grid
    openvdb::io::File(filename).write({grid});
}

/**
 * @brief      Load a binary file
 */
void ScalarField::load_binary() {
    std::ifstream infile(filename, std::ios::binary);

    // read data size
    this->gridsize = 1.0;
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            this->mat[i][j] = 0.0;
        }
        uint16_t nx = 0;
        infile.read((char*)&nx, sizeof(uint16_t));
        this->mat[i][i] = nx;
        this->grid_dimensions[i] = nx;
        this->gridsize *= nx;
    }

    this->calculate_inverse();
    this->calculate_volume();

    // prepare vector
    this->gridptr.resize(this->gridsize);

    // read float type
    uint16_t nx = 0;
    infile.read((char*)&nx, sizeof(uint16_t));
    if(nx == sizeof(float)) {   // directly load into vector
        infile.read((char*)&this->gridptr[0], sizeof(float) * this->gridsize);
    } else if(nx == sizeof(double)) {   // convert to double and load
        double val = 0.0;
        std::cout << "Reading " << this->gridsize << " values." << std::endl;
        for(unsigned int i=0; i<this->gridsize; i++) {
            infile.read((char*)&val, sizeof(double));
            this->gridptr[i] = (float)val;
        }
    } else {    // throw error
        throw std::runtime_error("Invalid data type when reading binary file.");
    }

    infile.close();

    this->scalar = 1.0;

    this->has_read = true;
    this->header_read = true;
}

/**
 * @brief      Load a cube file
 *
 * Technical details taken from:
 * http://paulbourke.net/dataformats/cube/
 */
void ScalarField::load_cube_file() {
    std::ifstream infile(filename);

    // convert bohr to angstrom
    static const float bohr_to_angstrom = 0.529177249;

    // string to read lines to
    std::string line;

    // skip the first two line
    for(unsigned int i=0; i<2; i++) {
        std::getline(infile, line);
    }

    //read number of atoms
    std::getline(infile, line);
    std::vector<std::string> pieces;
    boost::trim(line);
    unsigned nr_atoms = 0;
    try {
        boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
        nr_atoms = boost::lexical_cast<unsigned int>(pieces[0]);
        trans[0] = boost::lexical_cast<float>(pieces[1]) * bohr_to_angstrom;
        trans[1] = boost::lexical_cast<float>(pieces[2]) * bohr_to_angstrom;
        trans[2] = boost::lexical_cast<float>(pieces[3]) * bohr_to_angstrom;
    } catch(const std::exception& e) {
        std::cerr << "Cannot read number of atoms from " << pieces[0] << std::endl;
        std::cerr << "Line reads: \"" << line << "\"" << std::endl;
        exit(-1);
    }

    // grab size of box
    for(unsigned int i=0; i<3; i++) {
        // set matrix to zero
        for(unsigned int j=0; j<3; j++) {
            this->mat[i][j] = 0.0;
        }

        std::getline(infile, line);
        try {
            boost::trim(line);
            boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
            this->mat[i][i] = boost::lexical_cast<float>(pieces[0]) * boost::lexical_cast<float>(pieces[i+1]) * bohr_to_angstrom;
            this->grid_dimensions[i] = boost::lexical_cast<unsigned int>(pieces[0]);
            this->gridsize *= this->grid_dimensions[i];
        } catch(const std::exception& e) {
            std::cerr << "Cannot extract box size." << std::endl;
            std::cerr << "Line reads: \"" << line << "\"" << std::endl;
            exit(-1);
        }
    }

    this->calculate_inverse();
    this->calculate_volume();

    // prepare vector
    this->gridptr.reserve(this->gridsize);

    // skip atom lines
    for(unsigned int i=0; i<nr_atoms; i++) {
        std::getline(infile, line);
    }

    float_parser p;

    /* read spin up */
    unsigned int linecounter=0; // for the counter

    while(std::getline(infile, line)) {
        // set iterators
        std::string::const_iterator b = line.begin();
        std::string::const_iterator e = line.end();

        // parse
        std::vector<float> floats;
        boost::spirit::qi::phrase_parse(b, e, p, boost::spirit::ascii::space, floats);

        // expand gridptr with the new size
        unsigned int cursize = this->gridptr.size();
        this->gridptr.resize(cursize + floats.size());

        for(unsigned int j=0; j<floats.size(); j++) {
            this->gridptr[cursize + j] = floats[j];
        }

        linecounter++;
    }

    // Cube files are written in NX > NY > NZ ordering, whereas ScalarField requires NZ > NY > NX
    // thus we need to reorder the grid
    std::vector<float> newgrid(this->gridptr.size());
    for(unsigned int i=0; i<this->grid_dimensions[0]; i++) {    // x
        for(unsigned int j=0; j<this->grid_dimensions[1]; j++) {    // y
            for(unsigned int k=0; k<this->grid_dimensions[2]; k++) {    // z
                unsigned int idx_sf = k * this->grid_dimensions[0] * this->grid_dimensions[1] + j * this->grid_dimensions[0] + i;
                unsigned int idx_cub = i * this->grid_dimensions[1] * this->grid_dimensions[2] + j * this->grid_dimensions[2] + k;
                newgrid[idx_sf] = this->gridptr[idx_cub];
            }
        }
    }
    this->gridptr = newgrid;

    infile.close();

    this->scalar = 1.0;

    this->has_read = true;
    this->header_read = true;

    std::cout << "Done reading Gaussian Cube file" << std::endl;
    std::cout << "Read " << this->gridptr.size() << " values." << std::endl;
}
