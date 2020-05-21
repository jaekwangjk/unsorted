
#include<iostream> 
#include "heap_int.h"
// Driver program to test above functions 
int main() 
{ 
    MinHeap h(11);
    h.insertKey(3); 
    h.insertKey(2);  // value "key equals to value"
    //h.deleteKey(1); //index
    h.insertKey(15); 
    h.insertKey(5); 
    h.insertKey(4); 
    h.insertKey(45); 
    cout << h.extractMin() << " "; 
    cout << h.getMin() << " "; 
    h.decreaseKey(2, 1); 
    cout << h.getMin();
    
    
    int j;
    int m=30;
    int n=31;
    std::cout << std::endl; 
    j=m/2;
    std::cout << j <<" " ;
    j=n/2;
    std::cout << j <<" " ;
    
    return 0; 
} 
