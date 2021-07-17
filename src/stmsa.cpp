// stmsa.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

// #ifdef _DEBUG

// #define _CRTDBG_MAP_ALLOC
// #include <cstdlib>
// #include <crtdbg.h>

// the following two lines will make operator new tell crt where the leaked memory is allocated
// note that reloaded operator new cannot compile if the macros below are defined
// #define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
// #define new DEBUG_NEW

// #endif

#include "StarAlignment/StarAligner.hpp"
#include "Utils/Fasta.hpp"
#include "Utils/Utils.hpp"
#include "Utils/Arguments.hpp"

#include <boost/program_options.hpp>
#include <fstream>

void print_license();
bool parse_arguments(int argc, char **argv);

int main(int argc, char **argv)
{
#ifdef _DEBUG
    _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

    parse_arguments(argc, argv);

    auto start_point = std::chrono::system_clock::now();

    std::ifstream ifs(arguments::in_file_name);
    if (!ifs)
    {
        std::cout << "cannot access file " << arguments::in_file_name << '\n';
        exit(0);
    }

    const auto pseudo_sequences = utils::read_to_pseudo(ifs);
    std::cout << pseudo_sequences.size() << " sequences found\n";
    if (pseudo_sequences.size() < 2) exit(0);

    auto align_start = std::chrono::system_clock::now();
    auto insertions = star_alignment::StarAligner::get_gaps(pseudo_sequences);
    utils::print_duration(align_start, "aligning consumes"); std::cout << '\n';

    std::ofstream ofs(arguments::out_file_name);
    if (!ofs)
    {
        std::cout << "cannot write file " << arguments::out_file_name << '\n';
        exit(0);
    }

    ifs.clear();
    ifs.seekg(0);
    utils::insert_and_write(ofs, ifs, insertions);

    ifs.close();
    ofs.close();

    utils::print_duration(start_point, "total"); std::cout << '\n';
    print_license();

    return 0;
}

bool parse_arguments(int argc, char **argv)
{
    auto in_file_name_value = boost::program_options::value<std::string>(&arguments::in_file_name);
    in_file_name_value->required();

    auto out_file_name_value = boost::program_options::value<std::string>(&arguments::out_file_name);
    out_file_name_value->default_value("aligned.fasta");

    boost::program_options::options_description od("options");
    od.add_options()
        ("help,h", "produce help message")
        ("in,i", in_file_name_value, "name of the input file in fasta format")
        ("out,o", out_file_name_value, "name of the output file in fasta format")
        ("matrix", "output a character matrix rather than a fasta file");

    try
    {
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, od), vm);

        if (vm.find("help") != vm.cend()) { std::cout << od << '\n'; exit(0); }
        boost::program_options::notify(vm);

        arguments::output_matrix = vm.find("matrix") != vm.cend();
    }
    catch(const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        exit(0);
    }

    return true;
}

void print_license()
{
    std::cout << "http://lab.malab.cn/soft/halign/\n";
}
