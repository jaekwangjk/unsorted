#ifndef _Array_h_
#define _Array_h_

//////// product of array elements ////////
template<std::size_t Size>
auto mult(const std::array<double,Size>& a, const std::array<double,Size>& b) -> const std::array<double,Size>
{
    std::array<double,Size> c;
    for (int i=0;i<Size;i++)
    {
        c[i] = a[i]*b[i];
    }
    return c;
}
#endif

