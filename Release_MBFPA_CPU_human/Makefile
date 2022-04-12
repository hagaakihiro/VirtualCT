TARGET := mvct.exe
gC := g++
#CC := nvcc -arch=sm_20
#CC := nvcc -D_FORCE_INLINES -arch=sm_35
#CC := nvcc -arch=sm_35
#CC := nvcc -arch=sm_35 -DGPU -D_FORCE_INLINES -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -lm -lcudart -lcuda -lcublas -lcurand #-I/usr/local/cuda/include
op1 := -O2
op2 := -fopenmp


$(TARGET):		main.o	\
			mallocMD.o \
			reprojection_TOMO.o \
		     	lengthproj_TOMO.o \
			misc_TOMO.o \
			projectionDataTOMO.o \
			filterSinogramTOMO.o \
			prior_weight_production.o \
			IR_ImageTOMO.o
	$(gC)		main.o	\
			mallocMD.o \
		     	reprojection_TOMO.o \
		     	lengthproj_TOMO.o \
			misc_TOMO.o \
			projectionDataTOMO.o \
			filterSinogramTOMO.o \
			prior_weight_production.o \
			IR_ImageTOMO.o \
		$(op1) $(op2) \
		-o $(TARGET)
#		$(op1) $(op2) -L/usr/local/cuda/lib64 -lm -lcudart  -lcuda -lcublas -lcurand \
#		$(op1) $(op2)  -L/usr/local/cuda/lib64 -lm -lcudart\
#		     	projOnGPU_TOMO.o \
#		     	projOnGPU_TOMO.o \

main.o :		main.cpp \
			physParamsTOMO.h \
			mallocMD.h \
		     	lengthproj_TOMO.h \
		     	reprojection_TOMO.h \
			misc_TOMO.h \
			projectionDataTOMO.h \
			filterSinogramTOMO.h \
			virtual_projection.h \
			IR_ImageTOMO.h
	$(gC)	main.cpp \
		$(op1) $(op2) -c
#		     	projOnGPU_TOMO.h \

reprojection_TOMO.o:	reprojection_TOMO.cpp \
			reprojection_TOMO.h \
			physParamsTOMO.h \
			mallocMD.h
	$(gC)	reprojection_TOMO.cpp \
		$(op1) $(op2) -c

lengthproj_TOMO.o:	lengthproj_TOMO.cpp \
			lengthproj_TOMO.h \
			physParamsTOMO.h \
			mallocMD.h
	$(gC)	lengthproj_TOMO.cpp \
		$(op1) $(op2) -c

#projOnGPU_TOMO.o:	projOnGPU_TOMO.cu \
			projOnGPU_TOMO.h \
			physParamsTOMO.h \
			mallocMD.h
#	$(CC)	projOnGPU_TOMO.cu \
		$(op1) -c 

mallocMD.o:		mallocMD.cpp \
			mallocMD.h
	$(gC)	mallocMD.cpp \
		$(op1) $(op2) -c

misc_TOMO.o:		misc_TOMO.cpp \
			misc_TOMO.h \
			mallocMD.h
	$(gC)	misc_TOMO.cpp \
		$(op1) $(op2) -c

projectionDataTOMO.o:	projectionDataTOMO.cpp \
			projectionDataTOMO.h \
			physParamsTOMO.h \
			mallocMD.h
	$(gC)	projectionDataTOMO.cpp \
		$(op1) $(op2) -c

filterSinogramTOMO.o:	filterSinogramTOMO.cpp \
			filterSinogramTOMO.h \
			physParamsTOMO.h \
			mallocMD.h
	$(gC)	filterSinogramTOMO.cpp \
		$(op1) $(op2) -c

IR_ImageTOMO.o:	IR_ImageTOMO.cpp \
				IR_ImageTOMO.h \
				physParamsTOMO.h \
				mallocMD.h
	$(gC)	IR_ImageTOMO.cpp \
		$(op1) $(op2) -c

prior_weight_production.o: prior_weight_production.cpp \
		physParamsTOMO.h \
		mallocMD.h \
		virtual_projection.h
	$(gC)	prior_weight_production.cpp \
		$(op1) $(op2) -c
clean:
	rm -f *.o $(TARGET)
