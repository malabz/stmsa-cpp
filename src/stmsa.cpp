// stmsa.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#ifdef _DEBUG

#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

// the following two lines will make operator new tell crt where the leaked memory is allocated
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW

#endif

#include "UnitTest/PairwiseAlignmentTester.hpp"
#include "UnitTest/SuffixTreeTester.hpp"
#include "StarAlignment/StarAligner.hpp"
#include "Utils/Fasta.hpp"
#include "Utils/Utils.hpp"
#include "DAGLongestPath/DAGLongestPath.hpp"

#include <fstream>


int main(int argc, char **argv)
{
#ifdef _DEBUG
    _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif
    using namespace std::chrono;
 
    utils::Fasta fasta("D:\\data-set\\gisaid_hcov-19_2020_09_01_08_tmp.fasta");

    std::vector<std::vector<unsigned char>> pseudos; pseudos.reserve(fasta.sequences.size());
    utils::transform_to_pseudo(fasta.sequences.cbegin(), fasta.sequences.cend(), std::back_insert_iterator(pseudos));

    auto time_point = system_clock::now();
    auto aligned = star_alignment::StarAligner::align(pseudos);
    utils::print_duration(time_point, "aligning consumes"); std::cout << '\n';

    std::vector<std::string> outputs; outputs.reserve(fasta.sequences.size());
    utils::transform_from_pseudo(aligned.cbegin(), aligned.cend(), std::back_insert_iterator(outputs));

    std::ofstream ofs("D:\\tmp.fasta");
    utils::Fasta::write_to(ofs, outputs.cbegin(), outputs.cend());

    // std::vector<std::vector<unsigned>> graph{ { 0, 3, 2, 0, 0 }, { 0, 0, 0, 5, 0 }, { 0, 0, 0, 4, 8 }, { 0, 0, 0, 0, 1 }, { 0, 0, 0, 0, 0 } };
    // auto result = dag_longest_path::longest_path_of(graph, graph.size());
    // std::copy(result.cbegin(), result.cend(), std::ostream_iterator<size_t>(std::cout, ", "));
}
