/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   EDP is free software:                                                *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   EDP is distributed in the hope that it will be useful,               *
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
 * @param[in]  _filename        url to filename
 * @param[in]  _flag_is_locpot  whether this file is a locpot
 * @param[in]  _is_bin          is binary file
 */
ScalarField::ScalarField(const std::string &_filename, ScalarFieldInputFileType file_type) :
    filename(_filename),
    trans(0.0, 0.0, 0.0)
{
    if (!boost::filesystem::exists(this->filename)) {
        throw std::runtime_error("Cannot open " + this->filename + "!");
    }

    switch(file_type) {
        case ScalarFieldInputFileType::SFF_CHGCAR:
            this->scalar = -1;
            this->vasp5_input = false;
            this->has_read = false;
            this->header_read = false;
            this->flag_is_locpot = false;
        break;
        case ScalarFieldInputFileType::SFF_PARCHG:
            this->scalar = -1;
            this->vasp5_input = false;
            this->has_read = false;
            this->header_read = false;
            this->flag_is_locpot = false;
        break;
        case ScalarFieldInputFileType::SFF_LOCPOT:
            this->scalar = -1;
            this->vasp5_input = false;
            this->has_read = false;
            this->header_read = false;
            this->flag_is_locpot = true;
        break;
        case ScalarFieldInputFileType::SFF_CUB:
            this->load_cube_file();
        break;
        case ScalarFieldInputFileType::SFF_D2O:
            this->load_d2o_binary();
        break;
        default:
            throw std::runtime_error("Unknown file type.");
        break;
    }
}

#ifdef MOD_OPENVDB
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
                        accessor.setValue(ijk, std::abs(this->gridptr[idx]));
                    break;
                    case OpenVDB_METHOD::SQUARED:
                        accessor.setValue(ijk, this->gridptr[idx] * this->gridptr[idx]);
                    break;
                    case OpenVDB_METHOD::ABSOLUTE_LOG:
                        accessor.setValue(ijk, std::log(std::abs(this->gridptr[idx])));
                    break;
                    case OpenVDB_METHOD::POSITIVE:
                        accessor.setValue(ijk, std::max((fpt)0.0, this->gridptr[idx]));
                    break;
                    case OpenVDB_METHOD::NEGATIVE:
                        accessor.setValue(ijk, -std::min((fpt)0.0, this->gridptr[idx]));
                    break;
                    case OpenVDB_METHOD::POSITIVE_LOG:
                        accessor.setValue(ijk, std::log(std::max((fpt)1e-12, this->gridptr[idx])));
                    break;
                    case OpenVDB_METHOD::NEGATIVE_LOG:
                        accessor.setValue(ijk, -std::log(std::min((fpt)-1e-12, this->gridptr[idx])));
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
#endif // MOD_OPENVDB

/*
 * void output()
 *
 * Outputs a summary of the ScalarField to std::cout.
 * Mainly used for debugging purposes.
 *
 */
void ScalarField::output() const {
    std::cout << "Scalar: "<< this->scalar << std::endl;
    std::cout << std::endl;
    std::cout << "Matrix: ";
    for(unsigned i=0; i<3; i++) {
        for(unsigned j=0; j<3; j++) {
            std::cout << this->mat(i,j) << "\t";
        }
        std::cout << std::endl;
        std:: cout << "\t";
    }
    std::cout << std::endl;
    std::cout << "Inverse: ";
    for(unsigned i=0; i<3; i++) {
        for(unsigned j=0; j<3; j++) {
            std::cout << this->imat(i,j) << "\t";
        }
        std::cout << std::endl;
        std:: cout << "\t";
    }
    std::cout << std::endl;
    std::cout << "ion types: " << this->nrat.size() << " ( ";
    for(unsigned i=0; i<this->nrat.size(); i++) {
        std::cout << this->nrat[i] << " ";
    }
    std::cout << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "Grid dimensions: ";
    for(unsigned i=0; i<3; i++) {
        std::cout << this->grid_dimensions[i] << "\t";
    }
    std::cout << std::endl;
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

    // read scalar multiplication value
    try {
        this->read_scalar();
    } catch(const std::exception& e) {
        std::cout << "Error encountered in reading scalar value from CHGCAR" << std::endl;
        throw e;
    }

    // read matrix
    try {
        this->read_matrix();
    } catch(const std::exception& e) {
        std::cout << "Error encountered in reading unitcell matrix from CHGCAR" << std::endl;
        throw e;
    }

    // read number of atoms
    try {
        this->read_nr_atoms();
    } catch(const std::exception& e) {
        std::cout << "Error encountered in reading number of atoms from CHGCAR" << std::endl;
        throw e;
    }

    // read atomic positions
    try {
        this->read_atom_positions();
    } catch(const std::exception& e) {
        std::cout << "Error encountered in reading atomic positions from CHGCAR" << std::endl;
        std::cout << "ERROR: " << e.what() << std::endl;
        throw e;
    }

    // read grid dimensions
    try {
        this->read_grid_dimensions();
    } catch(const std::exception& e) {
        std::cout << "Error encountered in reading grid dimensions from CHGCAR" << std::endl;
        throw e;
    }

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
        this->scalar = boost::lexical_cast<fpt>(what[1]);
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
                mat(i,j) = boost::lexical_cast<fpt>(what[j+1]) * this->scalar;
            }
        }
    }

    // calculate matrix inverse and volume of unit cell
    this->imat = mat.inverse();
    this->volume = mat.determinant();
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
            boost::trim(line);
            // std::cout << line << std::endl;
            std::vector<std::string> pieces;
            boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
            // std::cout << pieces[0] << std::endl;
            // std::cout << pieces[1] << std::endl;
            // std::cout << pieces[2] << std::endl;

            fpt x = boost::lexical_cast<fpt>(pieces[0]);
            fpt y = boost::lexical_cast<fpt>(pieces[1]);
            fpt z = boost::lexical_cast<fpt>(pieces[2]);

            this->atom_pos.push_back(Vec3(x,y,z));
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
        std::vector<fpt> fpts;
        boost::spirit::qi::phrase_parse(b, e, p, boost::spirit::ascii::space, fpts);

        // expand gridptr with the new size
        unsigned int cursize = this->gridptr.size();
        this->gridptr.resize(cursize + fpts.size());

        // For CHGCAR type files, the electron density is multiplied by the cell volume
        // as described by the link below:
        // https://cms.mpi.univie.ac.at/vasp/vasp/CHGCAR_file.html#file-chgcar
        // Hence, for these files, we have to divide the value at the grid point by the
        // cell volume. For LOCPOT files, we should *not* do this.
        if(this->flag_is_locpot) {      // LOCPOT type files
            for(unsigned int j=0; j<fpts.size(); j++) {
                this->gridptr[cursize + j] = fpts[j];
            }
        } else {    // CHGCAR type files
            for(unsigned int j=0; j<fpts.size(); j++) {
                this->gridptr[cursize + j] = fpts[j] / this->volume;
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
 * fpt get_value_interp(x,y,z)
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
fpt ScalarField::get_value_interp(fpt x, fpt y, fpt z) const {
    if(!this->is_inside(x,y,z)) {
        return 0.0f;
    }

    // cast the input to grid space
    Vec3 r = this->realspace_to_grid(x,y,z);// - Vec3(0.5f, 0.5f, 0.5f);

    // recast
    if(r[0] < 0) r[0] += (fpt)this->grid_dimensions[0];
    if(r[1] < 0) r[1] += (fpt)this->grid_dimensions[1];
    if(r[2] < 0) r[2] += (fpt)this->grid_dimensions[2];

    // calculate value using trilinear interpolation
    fpt xd = remainderf(r[0], 1.0);
    fpt yd = remainderf(r[1], 1.0);
    fpt zd = remainderf(r[2], 1.0);

    if(xd < 0.0f) xd += 1.0f;
    if(yd < 0.0f) yd += 1.0f;
    if(zd < 0.0f) zd += 1.0f;

    fpt x0 = fmod(floor(r[0]), this->grid_dimensions[0]);
    fpt x1 = fmod(ceil(r[0]), this->grid_dimensions[0]);
    fpt y0 = fmod(floor(r[1]), this->grid_dimensions[1]);
    fpt y1 = fmod(ceil(r[1]), this->grid_dimensions[1]);
    fpt z0 = fmod(floor(r[2]), this->grid_dimensions[2]);
    fpt z1 = fmod(ceil(r[2]), this->grid_dimensions[2]);

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
bool ScalarField::is_inside(fpt x, fpt y, fpt z) const {
    // cast the input to the nearest grid point
    Vec3 d = this->realspace_to_direct(x,y,z);

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
 * fpt get_max_direction(dim)
 *
 * Get the maximum value in a particular dimension. This is a convenience
 * function for the get_value_interp() function.
 *
 */
fpt ScalarField::get_max_direction(unsigned int dim) {
    fpt sum = 0;
    for(unsigned int i=0; i<3; i++) {
        sum += this->mat(i,dim);
    }
    return sum;
}

/*
 * fpt get_value(i,j,k)
 *
 * Grabs the value at a particular grid point.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
fpt ScalarField::get_value(unsigned int i, unsigned int j, unsigned int k) const {
    unsigned int idx = k * this->grid_dimensions[0] * this->grid_dimensions[1] +
                       j * this->grid_dimensions[0] +
                       i;
    return this->gridptr[idx];
}

/*
 * Vec3 grid_to_realspace(i,j,k)
 *
 * Converts a grid point to a realspace vector. This function
 * is not being used at the moment.
 *
 */
Vec3 ScalarField::grid_to_realspace(fpt i, fpt j, fpt k) const {
    fpt dx = (fpt)i / (fpt)grid_dimensions[0];
    fpt dy = (fpt)j / (fpt)grid_dimensions[1];
    fpt dz = (fpt)k / (fpt)grid_dimensions[2];

    Vec3 r;
    r[0] = this->mat(0,0) * dx + this->mat(1,0) * dy + this->mat(2,0) * dz;
    r[1] = this->mat(0,1) * dx + this->mat(1,1) * dy + this->mat(2,1) * dz;
    r[2] = this->mat(0,2) * dx + this->mat(1,2) * dy + this->mat(2,2) * dz;

    return r;
}

/*
 * Vec3 realspace_to_grid(i,j,k)
 *
 * Convert 3d realspace vector to a position on the grid. Non-integer
 * values (i.e. fpting point) are given as the result.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
Vec3 ScalarField::realspace_to_grid(fpt i, fpt j, fpt k) const {
    Vec3 r = this->realspace_to_direct(i,j,k);

    r[0] *= fpt(this->grid_dimensions[0]);
    r[1] *= fpt(this->grid_dimensions[1]);
    r[2] *= fpt(this->grid_dimensions[2]);

    return r;
}

/*
 * Vec3 realspace_to_direct(i,j,k)
 *
 * Convert 3d realspace vector to direct position.
 *
 */
Vec3 ScalarField::realspace_to_direct(fpt i, fpt j, fpt k) const {
    Vec3 r;
    r[0] = this->imat(0,0) * i + this->imat(0,1) * j + this->imat(0,2) * k;
    r[1] = this->imat(1,0) * i + this->imat(1,1) * j + this->imat(1,2) * k;
    r[2] = this->imat(2,0) * i + this->imat(2,1) * j + this->imat(2,2) * k;

    return r;
}

fpt ScalarField::get_max() const {
    return *std::max_element(this->gridptr.begin(), this->gridptr.end());
}

fpt ScalarField::get_min() const {
    return *std::min_element(this->gridptr.begin(), this->gridptr.end());
}

Vec3 ScalarField::get_atom_position(unsigned int atid) const {
    if(atid < this->atom_pos.size()) {
        return this->mat * this->atom_pos[atid];
    } else {
        throw std::runtime_error("Requested atom id lies outside bounds");
    }
}

/**
 * @brief      Load a binary file
 */
void ScalarField::load_d2o_binary() {
    // read file format token
    std::ifstream infile(filename, std::ios::binary);
    char buf[3];
    infile.read(buf, 3);

    // verify token id
    if(std::string(buf,3) != "D2O") {
        throw std::runtime_error("File does not start with token D2O.");
    }

    // verify protocol id
    uint32_t protocol_id = 0;
    infile.read((char*)&protocol_id, sizeof(uint32_t));

    // check token id
    if(protocol_id == 0 || protocol_id > 2) {
        throw std::runtime_error("Invalid protocol id for d2o file: " + std::to_string(protocol_id));
    }

    // read matrix
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            float val = 0.0;
            infile.read((char*)&val, sizeof(float));
            this->mat(i,j) = val;
        }
    }

    // create inverse and calculate volume
    this->imat = this->mat.inverse();
    this->volume = this->mat.determinant();

    // read number of data points
    for(unsigned int i=0; i<3; i++) {
        uint32_t sz = 0;
        infile.read((char*)&sz, sizeof(uint32_t));
        this->grid_dimensions[i] = sz;
    }

    // read floating point size
    uint8_t fptsz = 0;
    infile.read((char*)&fptsz, sizeof(uint8_t));
    std::cout << "Recognizing floating point size: " << (int)fptsz << " bytes." << std::endl;

    if(sizeof(fpt) != fptsz) {
        throw std::runtime_error("Invalid floating point size. Incompatible with current grid pointer floating point type.");
    }

    // prepare vector
    this->gridsize = this->grid_dimensions[0] * grid_dimensions[1] * grid_dimensions[2];
    this->gridptr.resize(this->gridsize);

    // read data size
    uint64_t compdatasize = 0;
    infile.read((char*)&compdatasize, sizeof(uint64_t));
    std::cout << "Reading " << compdatasize << " bytes from file." << std::endl;

    // read data
    char* data = new char[compdatasize];
    infile.read(data, compdatasize);

    // decompress
    std::istringstream compressed(std::string(data, compdatasize));
    std::cout << "Building decompressor" << std::endl;
    boost::iostreams::filtering_istreambuf in;

    switch(protocol_id) {
        case 1:
            // GZIP compression
            in.push(boost::iostreams::gzip_decompressor());
        break;
        case 2:
            in.push(boost::iostreams::lzma_decompressor());
        break;
        default:
            throw std::runtime_error("Invalid protocol id for d2o file: " + std::to_string(protocol_id));
        break;
    }

    std::cout << "Decompressed data" << std::endl;
    in.push(compressed);

    std::ostringstream origin;
    boost::iostreams::copy(in, origin);
    const std::string gridptrdata = origin.str();
    memcpy(&this->gridptr[0], gridptrdata.data(), this->gridptr.size() * sizeof(fpt));

    this->scalar = 1.0;
    this->has_read = true;
    this->header_read = true;

    // clean up
    delete[] data;
    infile.close();

    std::cout << "Done reading D2O binary file" << std::endl;
    std::cout << "Read " << this->gridptr.size() << " values." << std::endl;
}

/**
 * @brief      Write to a binary file
 */
void ScalarField::write_d2o_binary(const std::string filename, uint32_t protocol_id) {
    // write file format token
    std::ofstream outfile(filename, std::ios::binary);
    char buf[] = "D2O";
    outfile.write(buf, 3);

    if(protocol_id == 0 || protocol_id > 2) {
        throw std::runtime_error("Invalid protocol id for d2o file: " + std::to_string(protocol_id));
    }

    // write protocol id
    outfile.write((char*)&protocol_id, sizeof(uint32_t));

    // write unit cell matrix
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            float val = this->mat(i,j);
            outfile.write((char*)&val, sizeof(float));
        }
    }

    // write number of data points
    for(unsigned int i=0; i<3; i++) {
        uint32_t sz = this->grid_dimensions[i];
        outfile.write((char*)&sz, sizeof(uint32_t));
    }

    // write floating point size
    uint8_t fptsz = sizeof(fpt);
    outfile.write((char*)&fptsz, sizeof(uint8_t));

    std::cout << "Floating point size determined at: " << (int)fptsz << " bytes" << std::endl;

    // compressing stream
    size_t gridptrsz = this->gridptr.size() * sizeof(fpt);
    char* data = new char[gridptrsz];
    memcpy(data, &this->gridptr[0], gridptrsz);
    std::istringstream origin(std::string(data, gridptrsz));
    boost::iostreams::filtering_istreambuf in;

    switch(protocol_id) {
        case 1:
            // GZIP compression
            in.push(
                boost::iostreams::gzip_compressor(
                    boost::iostreams::gzip_params(
                        boost::iostreams::gzip::best_compression
                    )
                )
            );
            std::cout << "Using GZIP compression." << std::endl;
        break;
        case 2:
            in.push(
                boost::iostreams::lzma_compressor()
            );
            std::cout << "Using LZMA compression." << std::endl;
        break;
        default:
            throw std::runtime_error("Invalid protocol id for d2o file: " + std::to_string(protocol_id));
        break;
    }

    in.push(origin);

    // store compression
    std::ostringstream compressed;
    boost::iostreams::copy(in, compressed);

    // output to file
    const std::string griddata = compressed.str();
    uint64_t sz = griddata.size();
    std::cout << "Compressed data to " << sz << " bytes ("
        << (float)sz / gridptrsz * 100 << " %)" << std::endl;

    // write data size and data
    outfile.write((char*)&sz, sizeof(uint64_t));
    outfile.write(griddata.data(), griddata.size());

    // clean up
    delete[] data;
    outfile.close();
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
            this->mat(i,j) = 0.0;
        }

        std::getline(infile, line);
        try {
            boost::trim(line);
            boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
            this->mat(i,i) = boost::lexical_cast<float>(pieces[0]) * boost::lexical_cast<float>(pieces[i+1]) * bohr_to_angstrom;
            this->grid_dimensions[i] = boost::lexical_cast<unsigned int>(pieces[0]);
            this->gridsize *= this->grid_dimensions[i];
        } catch(const std::exception& e) {
            std::cerr << "Cannot extract box size." << std::endl;
            std::cerr << "Line reads: \"" << line << "\"" << std::endl;
            exit(-1);
        }
    }

    // create inverse and calculate volume
    this->imat = this->mat.inverse();
    this->volume = this->mat.determinant();

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
