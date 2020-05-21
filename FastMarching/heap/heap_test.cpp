
#include<iostream> 
#include "heap.h"
#include<array>


// Should this be considered??
//#include"heap_bck.h"
//#include </Users/jaekwangkim/Program/eigen/Eigen/Dense>
//using Eigen::MatrixXd;


// Driver program to test above functions


//word , key & root ?

using namespace std;

int main() 
{
    // 2d heap data?
    const int dim = 2;
    typedef double data_t;

    std::array<int,dim> index;
    MinHeap<data_t,dim> h(11); 

    index = {1,2}; //what this 'index' is doing ??
    //I need to understand what is 2d heap data!
    
    HeapData<data_t,dim>  k(3.5,index);

    // along the these lines, index is fixed
    // and keep substitutitng different values
    h.insertKey(k); 
    k.set_val(3.1);
    h.insertKey(k); 
    k.set_val(2.9);
    h.insertKey(k); 
    k.set_val(15); 
    h.insertKey(k); 
    k.set_val(5);  
    h.insertKey(k); 
    k.set_val(4);  
    h.insertKey(k); 
    k.set_val(45); 
    h.insertKey(k);
    // it is not eleven elements, will this mean
    // The rest of h(11) is zero?
    // Then they are always the minumum!
    
    std::cout << h.extractMin().get_val() << std::endl;
    
    //extractMin does not have loop?
    
    for(const auto& element : h.getMin().get_bckPtr())
        std::cout << element ;
    std::cout << std::endl << h.getSize();
    // cout << h.getMin() << " ";
    // h.decreaseKey(2, 1); 
    // cout << h.getMin();
    return 0; 
} 

// what does "Decreases" value does?
