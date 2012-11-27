//reading from a text file
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(){
	int i;
	char *inname = "test.txt";
   ifstream infile (inname);
	cout << "Opened " << inname << " for reading." <<endl;
	while (infile >> i){
   	cout << i <<endl;
	}

  
   return 0;
}
