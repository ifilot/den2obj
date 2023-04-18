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

#include <chrono>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#include "scalar_field.h"
#include "config.h"
#include "isosurface.h"
#include "isosurface_mesh.h"

int main(int argc, char* argv[]) {
    try {
        TCLAP::CmdLine cmd("Converts VASP density file to wavefront object.", ' ', PROGRAM_VERSION);

        //**************************************
        // declare values to be parsed
        //**************************************

        // input filename
        TCLAP::ValueArg<std::string> arg_input_filename("i","input","Input file (i.e. CHGCAR)",true,"CHGCAR","filename");
        cmd.add(arg_input_filename);

        // output filename
        TCLAP::ValueArg<std::string> arg_output_filename("o","filename","Filename to output to",true,"test.png","string");
        cmd.add(arg_output_filename);

        // isovalue
        TCLAP::ValueArg<double> arg_isovalue("v", "isovalue", "Isovalue", false, 0.01, "double");
        cmd.add(arg_isovalue);

        // whether to center the wavefront object
        TCLAP::SwitchArg arg_c("c","center","center structure", cmd, false);

        // whether to transform the data carrier to another format, useful for converting CHGCAR
        // to binary format or to OpenVDB format
        TCLAP::SwitchArg arg_t("t","transform","Store in new format", cmd, false);

        cmd.parse(argc, argv);

        //**************************************
        // Inform user about execution
        //**************************************
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Executing "<< PROGRAM_NAME << " v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author:  Ivo Filot <i.a.w.filot@tue.nl>" << std::endl;
        std::cout << "Website: https://den2obj.imc-tue.nl" << std::endl;
        std::cout << "Github:  https://github.com/ifilot/den2obj" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;

        //**************************************
        // parsing values
        //**************************************
        std::string input_filename = arg_input_filename.getValue();
        std::string output_filename = arg_output_filename.getValue();

        // check whether only a transformation is being asked, if so, stop here
        std::unique_ptr<ScalarField> sf;
        if(input_filename.substr(input_filename.size()-4) == ".cub") {
            std::cout << "Opening " << input_filename << " as Gaussian cube file." << std::endl;
            sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_CUB);
        } else if(input_filename.substr(input_filename.size()-4) == ".d2o") {
            std::cout << "Opening " << input_filename << " as D2O binary file";
            sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_D2O);
        } else if(input_filename.substr(0,6) == "CHGCAR") {
            std::cout << "Opening " << input_filename << " as " << input_filename.substr(0,6) << std::endl;
            sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_CHGCAR);
        } else if(input_filename.substr(0,6) == "PARCHG") {
            std::cout << "Opening " << input_filename << " as " << input_filename.substr(0,6) << std::endl;
            sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_PARCHG);
        } else if(input_filename.substr(0,6) == "LOCPOT") {
            std::cout << "Opening " << input_filename << " as " << input_filename.substr(0,6) << std::endl;
            sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_LOCPOT);
        }

        // keep track of time
        auto start = std::chrono::system_clock::now();

        // check whether only a conversion is requested, if not, continue
        if(arg_t.getValue()) {
            if(output_filename.substr(input_filename.size()-4) == ".d2o") {
                sf->write_d2o_binary(output_filename);
            }
            #ifdef MOD_OPENVDB
            else if(output_filename.substr(input_filename.size()-6) == ".openvdb") {
                sf->write_to_vdb(output_filename, OpenVDB_METHOD::ABSOLUTE);
            }
            #endif
            else {
                std::runtime_error("Cannot interpret output file format. Please specify a valid extension.");
            }

        } else {
            double isovalue = arg_isovalue.getValue();
            std::cout << "Using isovalue: " << isovalue << std::endl;

            sf->read();
            std::cout << "Lowest value in scalar field: " << sf->get_min() << std::endl;
            std::cout << "Highest value in scalar field: " << sf->get_max() << std::endl;

            // construct isosurface generator object
            IsoSurface is(sf.get());
            is.marching_cubes(isovalue);

            // store path to extract filename
            boost::filesystem::path path(input_filename);

            // construct mesh storage object
            IsoSurfaceMesh ism(sf.get(), &is);
            ism.construct_mesh(arg_c.getValue());

            if(output_filename.substr(output_filename.size()-4) == ".obj") {
                ism.write_obj(output_filename, path.filename().string(), path.filename().string());
            } else if(output_filename.substr(output_filename.size()-4) == ".ply") {
                ism.write_ply(output_filename, path.filename().string(), path.filename().string());
            } else {
                std::runtime_error("Cannot interpret output file format. Please specify a valid extension.");
            }
        }

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << "Done in " << elapsed_seconds.count() << " seconds." << std::endl << std::endl;

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
