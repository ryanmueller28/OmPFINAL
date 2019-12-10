bin/name mpi.mpi: mpi-bs.cpp
	mpic++ -std=c++11 -Wall -O3 mpi-bs.cpp -o bin/mpi-bs.mpi
clean:
	rm bin/mpi-bs.mpi
