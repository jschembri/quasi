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


clean:
	rm -f quasi
	rm -f finiteVolume
	rm -f finiteDifference
