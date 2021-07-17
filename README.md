# Multiple Sequence Alignment with Suffix Tree

This project is aimed at improving a part of [HAlign](https://github.com/malabz/HAlign), which could align multiple nucleotide sequences fast and accurately.

## Usage

the command below prints the usage on windows

```
stmsa-x.x.x-win-x64.exe --help
```

## Change Log

 2021-07-17

  improve the output format

 2021-07-15

  further decrease the space requirement

 2021-04-25

  self-defined allocator for suffix tree, which reduce the construction time to about 1/3 comparing with ::operator new

 2021-03-29

  released as a crude but usable msa tool

 2021-03-11

  generalised suffixtree and needleman-wunsch algorithms implemented

## Build

 msvc

 clang

## Dependencies

 Boost

## License

  stmsa is licensed under the [MIT license](https://github.com/malabz/stmsa-cpp/blob/main/LICENSE).
