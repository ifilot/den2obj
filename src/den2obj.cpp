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

#include "config.h"
#include "scalar_field.h"
#include "generator.h"
#include "isosurface.h"
#include "isosurface_mesh.h"

int main(int argc, char* argv[]) {
    try {
        TCLAP::CmdLine cmd("Converts VASP density file to wavefront object.", ' ', PROGRAM_VERSION);

        //**************************************
        // declare values to be parsed
        //**************************************

        // input filename
        TCLAP::ValueArg<std::string> arg_input_filename("i","input","Input file (i.e. CHGCAR)",false,"CHGCAR","filename");
        cmd.add(arg_input_filename);

        // output filename
        TCLAP::ValueArg<std::string> arg_output_filename("o","filename","Filename to output to",true,"test.png","string");
        cmd.add(arg_output_filename);

        // isovalue
        TCLAP::ValueArg<double> arg_isovalue("v", "isovalue", "Isovalue", false, 0.01, "double");
        cmd.add(arg_isovalue);

        // whether to produce both positive as well as negative isovalue plots
        TCLAP::SwitchArg arg_d("d","dual","produce positive and negative versions of the isovalue", cmd, false);

        // whether to center the wavefront object
        TCLAP::SwitchArg arg_c("c","center","center structure", cmd, false);

        // whether to transform the data carrier to another format, useful for converting CHGCAR
        // to binary format or to OpenVDB format
        TCLAP::SwitchArg arg_t("t","transform","Store in new format", cmd, false);

        // output filename
        TCLAP::ValueArg<std::string> arg_generator("g","dataset","Dataset name",false,"","string");
        cmd.add(arg_generator);

        // select compression algorithm
        TCLAP::ValueArg<std::string> arg_algo("a","algo","Compression algorithm",false,"","string");
        cmd.add(arg_algo);

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
        std::string output_filename = arg_output_filename.getValue();

        // keep track of time
        auto start = std::chrono::system_clock::now();

        // verify algorithm choice
        D2OFormat::CompressionAlgo algo_id = D2OFormat::CompressionAlgo::AUTO; // default is auto-select
        if(!arg_generator.getValue().empty() || arg_t.getValue()) {
            if(!arg_algo.getValue().empty()) {
                auto got = D2OFormat::algos.find(arg_algo.getValue());
                if(got == D2OFormat::algos.end()) {
                    std::cerr << "Invalid choice for compression algorithm. Valid options are: " << std::endl;
                    for(const auto& i : D2OFormat::algos) {
                        std::cerr << " * " << i.first << std::endl;
                    }
                    throw std::runtime_error("Invalid choice for compression algorithm: " + arg_algo.getValue());
                } else {
                    algo_id = got->second;
                }
            }
        }

        // check if a generation is requested
        if(!arg_generator.getValue().empty()) {
            Generator gen;

            if(output_filename.substr(output_filename.size()-4) != ".d2o") {
                throw std::runtime_error("Invalid extension for dataset generation. Filename has to end in .d2o.");
            }

            std::cout << "Building grid using dataset: " << arg_generator.getValue() << std::endl;
            gen.build_dataset(arg_generator.getValue(), output_filename, algo_id);
        } else {
            const std::string input_filename = arg_input_filename.getValue();
            if(input_filename.empty()) {
                throw std::runtime_error("No input file is specified.");
            }

            boost::filesystem::path path_input_filename(input_filename);
            const std::string base_filename = path_input_filename.filename().string();

            // check whether only a transformation is being asked, if so, stop here
            std::unique_ptr<ScalarField> sf;
            if(base_filename.substr(base_filename.size()-4) == ".cub") {
                std::cout << "Opening " << input_filename << " as Gaussian cube file." << std::endl;
                sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_CUB);
            } else if(base_filename.substr(base_filename.size()-4) == ".d2o") {
                std::cout << "Opening " << input_filename << " as D2O binary file" << std::endl;
                sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_D2O);
            } else if(base_filename.substr(0,6) == "CHGCAR") {
                std::cout << "Opening " << input_filename << " as " << input_filename.substr(0,6) << std::endl;
                sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_CHGCAR);
            } else if(base_filename.substr(0,6) == "PARCHG") {
                std::cout << "Opening " << input_filename << " as " << input_filename.substr(0,6) << std::endl;
                sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_PARCHG);
            } else if(base_filename.substr(0,6) == "LOCPOT") {
                std::cout << "Opening " << input_filename << " as " << input_filename.substr(0,6) << std::endl;
                sf = std::make_unique<ScalarField>(input_filename, ScalarFieldInputFileType::SFF_LOCPOT);
            } else {
                throw std::runtime_error("Cannot interpret input file format. Please check the filename.");
            }

            // check whether only a conversion is requested, if not, continue
            if(arg_t.getValue()) {
                if(output_filename.substr(output_filename.size()-4) == ".d2o") {
                    sf->read();
                    sf->write_d2o_binary(output_filename, algo_id);
                    std::cout << "Writing as D2O binary file." << std::endl;
                }
                #ifdef MOD_OPENVDB
                else if(output_filename.substr(output_filename.size()-4) == ".vdb") {
                    sf->read();
                    sf->write_to_vdb(output_filename, OpenVDB_METHOD::ABSOLUTE);
                    std::cout << "Writing as OpenVDB file." << std::endl;
                }
                #endif
                else {
                    std::runtime_error("Cannot interpret output file format. Please specify a valid extension.");
                }

            } else {
                fpt isovalue = arg_isovalue.getValue();

                // automatically convert the isovalue to a positive number if d is set
                if(arg_d.getValue()) {
                    isovalue = std::abs(isovalue);
                }

                std::cout << "Using isovalue: " << isovalue << std::endl;

                sf->read();
                std::cout << "Lowest value in scalar field: " << sf->get_min() << std::endl;
                std::cout << "Highest value in scalar field: " << sf->get_max() << std::endl;

                // construct isosurface generator object
                IsoSurface is(sf.get());
                is.marching_cubes(isovalue);

                // store path to extract filename
                boost::filesystem::path path(output_filename);

                // construct mesh storage object
                IsoSurfaceMesh ism(sf.get(), &is);
                ism.construct_mesh(arg_c.getValue());

                // insert a suffix "pos" before the extension if dual conversion is set
                if(arg_d.getValue()) {
                    boost::filesystem::path p(output_filename);
                    auto ext = p.extension();
                    auto stem = p.stem();

                    // insert "_pos"
                    output_filename = stem.string() + "_pos" + ext.string();
                }

                // write isosurface to file; automatically capture file type from extension
                ism.write_to_file(output_filename, path.filename().string(), path.filename().string());

                // insert a suffix "pos" before the extension if dual conversion is set
                if(arg_d.getValue()) {
                    boost::filesystem::path p(output_filename);
                    auto ext = p.extension();
                    auto stem = p.stem();

                    // insert "_neg"
                    output_filename = stem.string().substr(0,stem.string().size()-4) + "_neg" + ext.string();

                    IsoSurface is_neg(sf.get());
                    is_neg.marching_cubes(-isovalue);

                    IsoSurfaceMesh ism_neg(sf.get(), &is_neg);
                    ism_neg.construct_mesh(arg_c.getValue());

                    // write isosurface to file; automatically capture file type from extension
                    ism_neg.write_to_file(output_filename, path.filename().string(), path.filename().string());
                }
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
