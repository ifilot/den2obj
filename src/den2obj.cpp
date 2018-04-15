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
        TCLAP::ValueArg<float> arg_isovalue("v","isovalue","Isovalue",true,0.1,"float");
        cmd.add(arg_isovalue);

        // whether to center the wavefront object
        TCLAP::SwitchArg arg_c("c","center","center structure", cmd, false);

        cmd.parse(argc, argv);

        //**************************************
        // Inform user about execution
        //**************************************
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Executing "<< PROGRAM_NAME << " v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <i.a.w.filot@tue.nl>" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;

        //**************************************
        // parsing values
        //**************************************
        std::string input_filename = arg_input_filename.getValue();
        std::string output_filename = arg_output_filename.getValue();
        float isovalue = arg_isovalue.getValue();

        std::cout << "Using isovalue: " << isovalue << std::endl;

        auto start = std::chrono::system_clock::now();

        ScalarField sf(input_filename, false);
        sf.read();

        IsoSurface is(&sf);
        is.marching_cubes(isovalue);

        // store path to extract filename
        boost::filesystem::path path(input_filename);

        IsoSurfaceMesh ism(&sf, &is);
        ism.construct_mesh(arg_c.getValue());
        ism.write_obj(output_filename, path.filename().string(), path.filename().string());

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << "You can directly import this file into blender using File > Import > Wavefront." << std::endl;
        std::cout << "NOTE: Recommended Blender import settings: Z-UP and Y-FORWARD." << std::endl;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << "Done in " << elapsed_seconds.count() << " seconds." << std::endl;

        return 0;

    }  catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
