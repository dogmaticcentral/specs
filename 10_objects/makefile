CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../specs
OPTM =  -O


VPATH =                       ../02_parsers:../04_residue_scoring:../06_main:\
        ../01_data_structures:../03_output:../05_cluster_scoring: ../07_utils


PLATFORM = $(shell uname)



LIBS =  -lm -lc
INCLUDE = -I ../00_include


DEFINES =  
SOURCES.F = 
SOURCES.f = 
SOURCES.c = 	specs.c  specs_aa_code.c  specs_area.c specs_adj_matrix.c specs_afa_out.c  \
        specs_alignment.c specs_clustering_out.c  specs_cluster_weight.c  \
	specs_chi2.c  specs_consensus_upgma.c specs_determinants.c \
	specs_entropy.c specs_fraction.c  specs_hybrid.c  specs_int_trace.c \
	specs_logger.c  specs_mapping.c   specs_node.c    specs_nj.c\
	specs_output_epi_spec.c  specs_overlap.c  specs_read_tree.c \
	specs_patch_almt.c specs_pheno.c \
	specs_read_cmd_file.c specs_read_alignment.c specs_read_dssp.c \
	specs_read_epi.c  specs_read_pdb.c specs_read_score.c \
	specs_score2coverage.c specs_score_out.c \
	specs_path_dets_out.c \
        specs_scoring.c specs_seq_pw_dist.c  specs_set_similarity.c  \
	specs_tree.c  specs_tree_out.c  specs_upgma.c  specs_utils.c \
        specs_valdar.c



AUX =  makefile  al.h specs_alignment.h specs_geometry.h specs_utils.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


