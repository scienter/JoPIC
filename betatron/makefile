EXEC = show
CC = mpicc
OBJS = main.o findparam.o parameterSetting.o shareData.o solveSpectrum.o
INCL = parameter.h constants.h
LIBS = -lm -lhdf5
$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
