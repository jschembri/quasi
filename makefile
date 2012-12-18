CXX=g++
CXXFLAGS=-Wall -g

quasi:
	$(CXX) $(CXXFLAGS) -c quasi.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o quasi quasi.o constants.o
finiteDifference:
	$(CXX) $(CXXFLAGS) -c finiteDifference.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o finiteDifference finiteDifference.o constants.o

finiteVolume:
	$(CXX) $(CXXFLAGS) -c finiteVolume.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o finiteVolume finiteVolume.o constants.o

analytic:
	$(CXX) $(CXXFLAGS) -c analytic.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o analytic analytic.o constants.o

analytic_subsonic:
	$(CXX) $(CXXFLAGS) -c analytic_subsonic.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o analytic_subsonic analytic_subsonic.o constants.o


subsonic:
	$(CXX) $(CXXFLAGS) -c subsonic.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o subsonic subsonic.o constants.o

flux:
	$(CXX) $(CXXFLAGS) -c flux.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o flux flux.o constants.o

flux:
	$(CXX) $(CXXFLAGS) -c subsonic.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o subsonic subsonic.o constants.o


newFlux:
	$(CXX) $(CXXFLAGS) -c newFlux.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o newFlux newFlux.o constants.o

machineZero:
	$(CXX) $(CXXFLAGS) -c machine_zero.cpp
	$(CXX) $(CXXFLAGS) -o machine_zero machine_zero.o

residual:
	$(CXX) $(CXXFLAGS) -c flux_with_residual.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o flux_with_residual flux_with_residual.o constants.o

Roe:
#	$(CXX) $(CXXFLAGS) -I /usr/local/include/eigen/ -c Roe.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -I /home/numericalMethods/eigen/ -c Roe.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o Roe Roe.o constants.o


Roe_supersonic:
#	$(CXX) $(CXXFLAGS) -I /usr/local/include/eigen/ -c Roe.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -I /home/numericalMethods/eigen/ -c Roe_supersonic.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o Roe_supersonic Roe_supersonic.o constants.o

pressure:
	$(CXX) $(CXXFLAGS) -c pressure.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o pressure pressure.o constants.o

pressure:
	$(CXX) $(CXXFLAGS) -c test.cpp 
	$(CXX) $(CXXFLAGS) -o test test.o 

clean:
	rm -f quasi
	rm -f finiteVolume
	rm -f finiteDifference
