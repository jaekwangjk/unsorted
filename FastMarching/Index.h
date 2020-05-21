#ifndef _Index_h_
#define _Index_h_

#include "Array.h"

template<typename=void>
constexpr int product()
{
    return 1;
}

template<typename... Ts>
constexpr int product(const int n, Ts... T)
{
    return n*product<>(T...);
}


// A template function to compute the indices from a flat index
//
// function template specialization - Does nothing
template<typename=void>
std::array<int,0> index(int K)
{
    return std::array<int,0>();
}

template<typename...T>
auto index(int K,int N, T...Ts) -> std::array<int,sizeof...(Ts) + 1>
{
    static constexpr std::size_t Nargs = sizeof...(Ts) + 1;
    std::array<int,Nargs>   local_ind;
    std::array<int,Nargs-1> local_partial_ind;
    local_ind[0] = K/product<>(Ts...);
    K = K - local_ind[0]*product<>(Ts...);

    local_partial_ind = index(K,Ts...);
    for (int i=1; i<Nargs; i++)
        local_ind[i] = local_partial_ind[i-1];
    return local_ind;
}

template<std::size_t Size, std::size_t...Indices>
std::array<int,Size> index_helper(int k,std::array<int,Size> arr,std::index_sequence<Indices...>)
{
    return index(k,std::get<Indices>(arr)...);
}

template<std::size_t Size>
std::array<int,Size> index(int k,std::array<int,Size> arr)
{
    return index_helper(k,arr,std::make_index_sequence<Size>());
}

#endif
