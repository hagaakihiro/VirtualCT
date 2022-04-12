TARGET := mvct.exe
CC := g++
#CC := nvcc
#CC := icc

$(TARGET):		main.o	\
			mallocMD.o \
			misc.o \
			projectionDataTOMO.o \
			filterSinogramTOMO.o \
                     reconstructImageTOMO.o
	$(CC)		main.o	\
			mallocMD.o \
			misc.o \
			projectionDataTOMO.o \
			filterSinogramTOMO.o \
			reconstructImageTOMO.o \
		-O2 -lm \
		-o $(TARGET)


main.o :		main.cpp \
			physParamsTOMO.h \
			mallocMD.h \
			misc.h \
			projectionDataTOMO.h \
			filterSinogramTOMO.h \
			reconstructImageTOMO.h
	$(CC)	main.cpp \
		-O2 -c

mallocMD.o:		mallocMD.cpp \
			mallocMD.h
	$(CC)	mallocMD.cpp \
		-O2 -c

misc.o:		misc.cpp \
			misc.h \
			mallocMD.h
	$(CC)	misc.cpp \
		-O2 -c

projectionDataTOMO.o:	projectionDataTOMO.cpp \
			projectionDataTOMO.h \
			physParamsTOMO.h \
			mallocMD.h
	$(CC)	projectionDataTOMO.cpp \
		-O2 -c

filterSinogramTOMO.o:	filterSinogramTOMO.cpp \
			filterSinogramTOMO.h \
			physParamsTOMO.h \
			mallocMD.h
	$(CC)	filterSinogramTOMO.cpp \
		-O2 -c

reconstructImageTOMO.o:	reconstructImageTOMO.cpp \
				reconstructImageTOMO.h \
				physParamsTOMO.h \
				mallocMD.h
	$(CC)	reconstructImageTOMO.cpp \
		-O2 -c

clean:
	rm -f *.o $(TARGET)
