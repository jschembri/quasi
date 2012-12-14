#include <iostream>

using namespace std;


int * fillarr(int arr[], int length){
   for (int i = 0; i < length; ++i){
   	arr[i] = 1;
   }
   return arr;
}


int main(int argc, char **argv){
	int arr[5];
	for (int i=0;i<=4;i++){
		arr[i] = 0;
	}
	fillarr(arr,5);	

	for (int i=0;i<=4;i++){
		cout << arr[i] << endl;
	}

	return 0;
}


