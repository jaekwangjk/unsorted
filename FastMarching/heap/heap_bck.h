// A C++ program to demonstrate common Binary Heap Operations 
// Need to check how this heap data structure scales
#ifndef HEAP_H
#define HEAP_H

#include<iostream> 
#include<climits> 
#include<array> 
#include<vector> 
  
// A class for Min Heap 

// Ideal for template class with template parameters: dim,

template<class T, int dim>
class HeapData
{
private:
    T val;
    std::array<int,dim> bckPtr;
public: 
    // constructor
    HeapData(T v, std::array<int,dim> ptr):val(v),bckPtr(ptr) { }

    HeapData(): val(),bckPtr() {}

    T get_val() {return val;}

    std::array<int,dim>& get_bckPtr() {return this->bckPtr;}

    void set_val(T v) {val = v;}

    void set_bckPtr(std::array<int,dim> ind) {bckPtr = ind;}

    friend bool operator < (const HeapData<T,dim>& data1, const HeapData<T,dim>& data2)
    {
        bool res;
        return (data1.val < data2.val)? true : false;

    }
    friend bool operator > (const HeapData<T,dim>& data1, const HeapData<T,dim>& data2)
    {
        bool res;
        //return (data1.val > data2.val)? true : false;
        return (data1.val > data2.val)? true : false;
    }
};


template<class T,int dim>
struct MinHeap 
{ 
    std::vector<HeapData<T,dim>> harr;  // Array of elements in heap 
    Matrix<dim,int>:: locations;         // Array to describe the position 
    int  heap_size;                     // Current number of elements in min heap 

    // Constructor: Builds a heap from a given array a[] of given size 
    MinHeap(const int& nx, const int& ny):heap_size(0),harr(nx*ny),locations(nx,ny){ } 
  
    // A recursive method to heapify a subtree with the root at given index 
    // This method assumes that the subtrees are already heapified 
    //
    void MinHeapify(int i) 
    { 
        int l = left(i); 
        int r = right(i); 
        int smallest = i; 
        if (l < heap_size && harr[l] < harr[i]) 
            smallest = l; 
        if (r < heap_size && harr[r] < harr[smallest]) 
            smallest = r; 
        if (smallest != i) 
        { 
            std::swap(harr[i], harr[smallest]); 
            locations(harr.
            MinHeapify(smallest); 
        } 
    } 
  
    int parent(int i) { return (i-1)/2; } 
  
    // to get index of left child of node at index i 
    int left(int i) { return (2*i + 1); } 
  
    // to get index of right child of node at index i 
    int right(int i) { return (2*i + 2); } 
  
    // Method to remove minimum element (or root) from min heap 
    HeapData<T,dim> extractMin() 
    { 
        if (heap_size <= 0) 
        {
            HeapData<T,dim> root;
            root.set_val(INT_MAX);
            return root;
        }
        if (heap_size == 1) 
        { 
            heap_size--; 
            return harr[0]; 
        } 
      
        // Store the minimum value, and remove it from heap 
        HeapData<T,dim> root = harr[0]; 
        harr[0] = harr[heap_size-1]; 
        heap_size--; 
        MinHeapify(0); 
      
        return root; 
    } 
  
    // Decreases value of key at index 'i' to new_val.  It is assumed that 
    // new_val is smaller than harr[i]. 
    void decreaseKey(int i, T new_val) 
    { 
        harr[i].set_val(new_val); 
        while (i != 0 && harr[parent(i)] > harr[i]) 
        { 
            std::swap(harr[i], harr[parent(i)]); 
           i = parent(i); 
        } 
    } 
  
    // Returns the minimum key (key at root) from min heap 
    HeapData<T,dim> getMin() { return harr[0]; } 
  
    // This function deletes key at index i. It first reduced value to minus 
    // infinite, then calls extractMin() 
    void deleteKey(int i) 
    { 
        decreaseKey(i, INT_MIN); 
        extractMin(); 
    } 
  
    // Inserts a new key 'k' 
    void insertKey(HeapData<T,dim> k) 
    { 
        if (heap_size == harr.size()) 
        { 
            std::cout << "\nOverflow: Could not insertKey\n"; 
            return; 
        } 
      
        // First insert the new key at the end 
        heap_size++; 
        int i = heap_size - 1; 
        harr[i] = k; 
      
        // Fix the min heap property if it is violated 
        while (i != 0 && harr[parent(i)] > harr[i]) 
        { 
            std::swap(harr[i], harr[parent(i)]); 
            i = parent(i); 
        } 
    } 

}; 
  
  
#endif  
