FC = gfortran
FCFLAGS = -O3

OBJS = mod_math.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: Ionboost

Ionboost: IonBoost28.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

clean:
	$(RM) IonBoost28 *.o *.mod
