# Executable
EXE = program.exe

# Optimization
# OPT = -Ofast
OPT = -O2

# Compiler
CC = gcc
FC = gfortran

# Compiler Library
CLIB = -lm -lgfortran -lquadmath

# Compiler Flags
CFLG = -Wall -Wextra -pedantic $(OPT) -march=native
F90FLG = -Wall -Wextra -Ofast -static -march=native
F77FLG = -Ofast -static -march=native

# Directories
BDIR = bin
SDIR = src
HDIR = include
ODIR = .obj
DDIR = .dep

# Files
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
CDEP = $(patsubst $(SDIR)/%.c, $(DDIR)/%.d, $(CSRC))
F90SRC = $(wildcard $(SDIR)/*.f90)
F90OBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(F90SRC))
F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))

# Dependencies Flags
DFLG = -MMD -MF $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $@)

# Targets
all: $(BDIR) $(ODIR) $(DDIR) $(BDIR)/$(EXE)

$(BDIR) $(ODIR) $(DDIR):
	mkdir -p $@

$(BDIR)/$(EXE): $(COBJ) $(F90OBJ) $(F77OBJ)
	$(CC) -o $@ $^ -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLG) $(DFLG) -c $< -o $@ -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(F90FLG) -c $^ -o $@

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(F77FLG) -c $^ -o $@

# Clean
.PHONY: clean cleanall 
clean:
	$(RM) $(BDIR)/$(EXE) $(ODIR)/*.o $(DDIR)/*.d

cleanall: clean
	$(RM) -r $(BDIR) $(ODIR) $(DDIR)

FigureRefl:
	$(MAKE) -C Data/Reflection/
	xdg-open Data/Reflection/Figure.pdf

FigureTLGF:
	$(MAKE) -C Data/TLGF/

FigureSI:
	$(MAKE) -C Data/SI/

FigureChew:
	$(MAKE) -C Data/Chew/

FigureNadson:
	$(MAKE) -C Data/Nadson/

FigurePlaneWave:
	$(MAKE) -C Data/PlaneWave/

FigureFEKOPlaneWave:
	$(MAKE) -C Data/FEKOPlaneWave/

FigureFarFiled:
	$(MAKE) -C Data/FarField/

FigurePaulus:
	$(MAKE) -C Data/Paulus/

FigureBragg:
	$(MAKE) -C Data/Bragg/

-include $(CDEP)