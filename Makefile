CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
	exist = $(shell if [ -e '/usr/local/bin/g++-10' ]; then echo "exist"; else echo "notexist"; fi;)
	ifeq ($(exist),exist)
		CC:=g++-10
	else
        	exist = $(shell if [ -e '/usr/local/bin/g++-9' ]; then echo "exist"; else echo "notexist"; fi;)
        	ifeq ($(exist),exist)
                	CC:=g++-9
		else
			CC:=g++-8
		endif
	endif
endif

HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse
EXE_ICD=bin/RamEx-IRCA-clean
EXE_ICM=bin/RamEx-IRCA-cmt
EXE_IGC=bin/RamEx-IRCA-gcm
EXE_IGR=bin/RamEx-IRCA-group
EXE_IPC=bin/RamEx-IRCA-PCA
EXE_ISR=bin/RamEx-IRCA-srt
EXE_IWN=bin/RamEx-IRCA-wnt
EXE_RQC=bin/RamEx-QC
EXE_RRF=bin/RamEx-RBCS-RF

tax:
	$(CC) -o $(EXE_ICD) src/Ramex_IRCA_clean_data.cpp $(HASHFLG)	
	$(CC) -o $(EXE_ICM) src/Ramex_IRCA_correlation_matrix_transfer.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_IGC) src/Ramex_IRCA_Generate_Correlation_Matrix.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_IGR) src/Ramex_IRCA_Grouping.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_IPC) src/Ramex_IRCA_PCA.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_ISR) src/Ramex_IRCA_spec_range_trimming.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_IWN) src/Ramex_IRCA_wave_number_trimming.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_RQC) src/Ramex_QC.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_RRF) src/Ramex_RBCS_Peak.cpp $(HASHFLG) $(BUILDFLG)


clean:
	rm -rf bin/* src/*.o

