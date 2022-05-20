#include <iostream>

using namespace std;

int main(){

	int *A;
	A = new int[4]{0,1,4,5};

	for (int i =0; i<4;i++) cout<<A[i]<<" ";
	return 0;
}
