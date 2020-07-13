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
        TCLAP::ValueArg<float> arg_isovalue("v", "isovalue", "Isovalue",false, 0.01, "float");
        cmd.add(arg_isovalue);

        // whether to center the wavefront object
        TCLAP::SwitchArg arg_c("c","center","center structure", cmd, false);

        // whether input file is binary
        TCLAP::SwitchArg arg_b("b","binary","binary file", cmd, false);

        // whether to write to vdb file
        TCLAP::SwitchArg arg_d("d","vdb","vdb file", cmd, false);

        // whether to write to ply or to wavefront file
        TCLAP::SwitchArg arg_p("p","ply","ply file", cmd, false);

        // which openVDB method to use
        TCLAP::ValueArg<std::string> arg_method("m","method","Which OpenVDB method to use",false,"absolute","string");
        cmd.add(arg_method);

        cmd.parse(argc, argv);

        //**************************************
        // Inform user about execution
        //**************************************
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Executing "<< PROGRAM_NAME << " v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <i.a.w.filot@tue.nl>" << std::endl;
        std::cout << "Website: https://github.com/ifilot/den2obj" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;

        //**************************************
        // parsing values
        //**************************************
        std::string input_filename = arg_input_filename.getValue();
        std::string output_filename = arg_output_filename.getValue();
        float isovalue = arg_isovalue.getValue();

        std::cout << "Using isovalue: " << isovalue << std::endl;

        auto start = std::chrono::system_clock::now();

        if(arg_b.getValue()) {
            std::cout << "Opening " << input_filename << " as binary file." << std::endl;
        } else {
            if(input_filename.substr(input_filename.size()-4) == ".cub") {
                std::cout << "Opening " << input_filename << " as Gaussian cube file." << std::endl;
            } else {
                std::cout << "Opening " << input_filename << " as VASP CHGCAR file." << std::endl;
            }
        }

        ScalarField sf(input_filename, false, arg_b.getValue());
        sf.read();

        std::cout << "Lowest value in scalar field: " << sf.get_min() << std::endl;
        std::cout << "Highest value in scalar field: " << sf.get_max() << std::endl;

        if(arg_d.getValue()) {
            std::cout << "Creating OpenVDB file" << std::endl;
            std::cout << "Using method flag: " << arg_method.getValue() << std::endl;

            if(arg_method.getValue() == "absolute") {
                sf.write_to_vdb(output_filename, OpenVDB_METHOD::ABSOLUTE);
            } else if(arg_method.getValue() == "absolute_log") {
                sf.write_to_vdb(output_filename, OpenVDB_METHOD::ABSOLUTE_LOG);
            }  else if(arg_method.getValue() == "positive") {
                sf.write_to_vdb(output_filename, OpenVDB_METHOD::POSITIVE);
            }  else if(arg_method.getValue() == "negative") {
                sf.write_to_vdb(output_filename, OpenVDB_METHOD::NEGATIVE);
            }  else if(arg_method.getValue() == "positive_log") {
                sf.write_to_vdb(output_filename, OpenVDB_METHOD::POSITIVE_LOG);
            }  else if(arg_method.getValue() == "negative_log") {
                sf.write_to_vdb(output_filename, OpenVDB_METHOD::NEGATIVE_LOG);
            } else {
                throw std::runtime_error("Invalid method supplied: " + arg_method.getValue());
            }

            std::cout << "Output has been written to: " << arg_output_filename.getValue() << std::endl;

            return 0;
        }

        IsoSurface is(&sf);
        is.marching_cubes(isovalue);

        // store path to extract filename
        boost::filesystem::path path(input_filename);

        IsoSurfaceMesh ism(&sf, &is);
        ism.construct_mesh(arg_c.getValue());

        if(arg_p.getValue()) {
            ism.write_ply(output_filename, path.filename().string(), path.filename().string());
        } else {
            ism.write_obj(output_filename, path.filename().string(), path.filename().string());
        }

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        if(arg_p.getValue()) {
            std::cout << "--------------------------------------------------------------------------------------" << std::endl;
            std::cout << "You can directly import this file into blender using File > Import > Wavefront (.obj)." << std::endl;
            std::cout << "--------------------------------------------------------------------------------------" << std::endl;
        } else {
            std::cout << "-------------------------------------------------------------------------------------" << std::endl;
            std::cout << "You can directly import this file into blender using File > Import > Stanford (.ply)." << std::endl;
            std::cout << "NOTE: Recommended Blender import settings: Z-UP and Y-FORWARD." << std::endl;
            std::cout << "-------------------------------------------------------------------------------------" << std::endl;
        }
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << "Done in " << elapsed_seconds.count() << " seconds." << std::endl << std::endl;

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
