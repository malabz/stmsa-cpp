CC=			clang++
CFLAGS=		-g -Wall -Wc++-compat -std=c++17 -O3 -static
CPPFLAGS=	
INCLUDES=
OBJS=		Utils/Insertion.o Utils/Insertion.o Utils/Graph.o Utils/Graph.o Utils/Fasta.o Utils/Fasta.o Utils/Arguments.o Utils/Arguments.o Utils/Utils.o Utils/Utils.o Utils/NucleicAcidColumn.o Utils/NucleicAcidColumn.o \
		SuffixTree/LeftChildRightSiblingBinaryTree.o SuffixTree/LeftChildRightSiblingBinaryTree.o \
		PairwiseAlignment/PairwiseAlignment.o PairwiseAlignment/PairwiseAlignment.o \
		StarAlignment/StarAligner.o StarAlignment/StarAligner.o
PROG=		stmsa
LIBS=		-lboost_program_options

.SUFFIXES:.c .o
.PHONY:all clean


.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)


stmsa:$(OBJS) stmsa.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)


clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM $(OBJS)


