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

read:
	$(CXX) $(CXXFLAGS) -c read.cpp
	$(CXX) $(CXXFLAGS) -o read read.o

pressure:
	$(CXX) $(CXXFLAGS) -c pressure.cpp constants.cpp
	$(CXX) $(CXXFLAGS) -o pressure pressure.o constants.o

read2:
	$(CXX) $(CXXFLAGS) -c read2.cpp
	$(CXX) $(CXXFLAGS) -o read2 read2.o

clean:
	rm -f quasi
	rm -f finiteVolume
	rm -f finiteDifference
