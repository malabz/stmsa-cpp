# Multiple Sequence Alignment with Suffix Tree

This project is aimed at improving a part of [HAlign2.0](https://github.com/malabz/HAlign), which could align multiple nucleotide sequences fast and accurately.

## Usage

```
stmsa-x.x.x-win-x64.exe unaligned.fasta output.fasta
```

## Change Log

* 2021-04-25

  self-defined allocator for suffix tree, which reduce the construction time to about 1 / 3 comparing with ::operator new

* 2021-03-29

  released as a crude but usable msa tool

* 2021-03-11

  generalised suffixtree and needleman-wunsch algorithms implemented

## Build

- msvc

- clang

## License

  stmsa is licensed under the [MIT license](https://github.com/malabz/stmsa-cpp/blob/main/LICENSE).
