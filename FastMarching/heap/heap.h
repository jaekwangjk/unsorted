// A C++ program to demonstrate common Binary Heap Operations 
// Need to check how this heap data structure scales
#ifndef HEAP_H
#define HEAP_H

#include<iostream> 
#include<climits> 
#include<array> 
#include<vector> 
  
// A class for Min Heap 

template<class T> // This is standard approach for function template
class HeapData
{
protected:
    T val;  ///?? T::  why not "double val" or data_t val"
    int bckPtr;// back to grid!
public: 
    // constructor
    HeapData(T v, int ptr):val(v),bckPtr(ptr) { }

    HeapData(): val(),bckPtr() {}

    const T& get_val() {return val;}

    const int& get_bckPtr() {return this->bckPtr;}

    void set_val(T v) {val = v;}

    void set_bckPtr(int ind) {bckPtr = ind;}

    // Friend class allows "assess" to protected member...
    friend bool operator < (const HeapData<T>& data1, const HeapData<T>& data2)
    {
        bool res;
        return (data1.val < data2.val)? true : false;

    }
    friend bool operator > (const HeapData<T>& data1, const HeapData<T>& data2)
    {
        bool res;
        //return (data1.val > data2.val)? true : false;
        return (data1.val > data2.val)? true : false;
    }
};


template<class T>
struct MinHeap 
{
    // why class ? and struct?
    
    //(???) how this type would looks like?
    std::vector<HeapData<T>> harr;  // Array of elements in heap 
    //(????????) is harr size is already defined? or should be defined later???
    
    int  heap_size;                 // Current number of elements in min heap
    
    //(???) how does index and location different?
    std::vector<int> locations;

    // Constructor: Builds a heap from a given array a[] of given size 
    MinHeap(const int& pcount):heap_size(0),harr(pcount),locations(pcount,-1){}
    ///(???), what this constructor does this? why heap_size is intialized as 0?
    // is location == bakPrt? what pcount,-1?
    //locatin is kinda reverse of backptr
    
    // A recursive method to heapify a subtree with the root at given index 
    // This method assumes that the subtrees are already heapified 
    //
    // It seems MinHeatfy find min out of 3 (1-parents 2-children)
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
            std::swap(locations[harr[i].get_bckPtr()], locations[harr[smallest].get_bckPtr()]); 
            MinHeapify(smallest); 
        } 
    } 
  
    int parent(int i) { return (i-1)/2; } 
  
    // to get index of left child of node at index i 
    int left(int i) { return (2*i + 1); } 
  
    // to get index of right child of node at index i 
    int right(int i) { return (2*i + 2); } 
  
    // Method to remove minimum element (or root) from min heap 
    HeapData<T> extractMin() 
    { 
        if (heap_size <= 0) 
        {
            HeapData<T> root;
            root.set_val(INT_MAX);
            return root;
        }
        if (heap_size == 1) 
        { 
            heap_size--;
            //why (-1) here? location and bckPtr. how different?
            locations[harr[0].get_bckPtr()] = -1;
            return harr[0]; 
        } 
      
        // Store the minimum value, and remove it from heap 
        HeapData<T> root = harr[0]; 
        locations[root.get_bckPtr()] = -1; // help? what is doing this?

        harr[0] = harr[heap_size-1]; //what?? harr[0] root! but harr[heap_size-1] is
        // long after root!
        locations[harr[0].get_bckPtr()] = 0;
        heap_size--; 
        MinHeapify(0); //heapify the first zero root!
      
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
            std::swap(locations[harr[i].get_bckPtr()], locations(harr[parent(i)].get_bckPtr())); 
            i = parent(i); 
        } 
    } 
  
    // Returns the minimum key (key at root) from min heap 
    HeapData<T> getMin() { return harr[0]; } 
  
    // This function deletes key at index i. It first reduced value to minus 
    // infinite, then calls extractMin() 
    void deleteKey(int i) 
    { 
        decreaseKey(i, INT_MIN); 
        extractMin(); 
    } 
  
    // Inserts a new key 'k' 
    void insertKey(HeapData<T> k) 
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
        locations[k.get_bckPtr()] = i;
      
        // Fix the min heap property if it is violated 
        while (i != 0 && harr[parent(i)] > harr[i]) 
        { 
            std::swap(harr[i], harr[parent(i)]); 
            std::swap(locations[harr[i].get_bckPtr()], locations[harr[parent(i)].get_bckPtr()]); 
            i = parent(i); 
        } 
    } 

}; 
  
  
#endif  
