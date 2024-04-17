#include <den2obj/scalar_field.h>
#include <den2obj/generator.h>
#include <memory>
#include <string>

int main() {
    // construct scalar field
    Generator gen;
    const std::string filename = "genus2.d2o";
    gen.build_dataset("genus2", filename, D2OFormat::CompressionAlgo::BZIP2);

    auto sf = std::make_unique<ScalarField>(filename, 
                                            ScalarFieldInputFileType::SFF_D2O);

    return 0;
}