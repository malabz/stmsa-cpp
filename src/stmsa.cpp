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

#include <fstream>

void print_usage(const char* argv0);
void print_license();

int main(int argc, char **argv)
{
#ifdef _DEBUG
    _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

    if (argc != 3) { print_usage(argv[0]); exit(0); }

    auto start_point = std::chrono::system_clock::now();

    std::ifstream ifs(argv[1]);
    if (!ifs)
    {
        std::cout << "cannot access file " << argv[1] << '\n';
        exit(0);
    }

    const auto [identifications, sequences] = utils::read_to_pseudo(ifs);
    std::cout << sequences.size() << " sequences found\n";
    if (sequences.size() < 2) exit(0);

    auto align_start = std::chrono::system_clock::now();
    auto aligned = star_alignment::StarAligner::align(sequences);
    utils::print_duration(align_start, "aligning consumes"); std::cout << '\n';

    std::vector<std::string> outputs; outputs.reserve(sequences.size());
    utils::transform_from_pseudo(aligned.cbegin(), aligned.cend(), std::back_insert_iterator(outputs));

    std::ofstream ofs(argv[2]);
    if (!ofs)
    {
        std::cout << "cannot write file " << argv[2] << '\n';
        exit(0);
    }

    utils::Fasta::write_to(ofs, outputs.cbegin(), outputs.cend());

    utils::print_duration(start_point, "total"); std::cout << '\n';
    print_license();

    return 0;
}

void print_usage(const char* argv0)
{
    std::cout << "usage: \n";
    std::cout << argv0 << " unaligned.fasta output.fasta\n";
    print_license();
}

void print_license()
{
    std::cout << "http://lab.malab.cn/soft/halign/\n";
}
