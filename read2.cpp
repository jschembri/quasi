#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    double i;
	 int j=0;
    char *inname = "test.txt";
    ifstream infile(inname);
	 double pressure[5];

    if (!infile) {
        cout << "There was a problem opening file "
             << inname
             << " for reading."
             << endl;
        return 0;
    }
    cout << "Opened " << inname << " for reading." << endl;
    while (infile >> i) {
		  pressure[j] = i;
        cout << "Value from file is " << i << endl;
	     j +=1;
    }

	 for (int i=0;i<5;i++){
       cout << "The pressure array is: " << pressure[i] << endl;
	 }
    return 0;
}

