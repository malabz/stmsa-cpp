#include "Pseudo.h"

std::string pseudo::to_pseudo(const std::string& str)
{
    std::string pseu;
    for (auto i : str) pseu.push_back(to_pseudo(i));
    return pseu;
}

inline char pseudo::to_pseudo(char ch)
{
    switch (ch)
    {
    case '-':
        return GAP;
    case 'C':
    case 'c':
        return C;
    case 'G':
    case 'g':
        return G;
    case 'A':
    case 'a':
        return A;
    case 'T':
    case 't':
        return T;
    default:
        return UNKNOWN;
    }
}

