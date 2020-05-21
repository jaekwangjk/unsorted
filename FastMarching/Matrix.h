#ifndef _Matrix_h_
#define _Matrix_h_

#include <vector>       /* pow */

template <const int dim, typename T=double>
class Matrix : public std::vector<Matrix<dim-1,T>>
{
    public:
        template<typename...Ints>
        Matrix(const int& n,const Ints&...ints) : std::vector<Matrix<dim-1,T>>(n,Matrix<dim-1,T>(ints...))
        {
        }

        template<typename...Ints>
        const T& operator ()(const int& i, const Ints&...ints) const
        {
            return this->operator[](i)(ints...);
        }

        template<typename...Ints>
        T& operator ()(const int& i, const Ints&...ints)
        {
            return this->operator[](i)(ints...);
        }

};

template<typename T>
class Matrix<1,T> : public std::vector<T>
{
    public:
        Matrix(const int& n) : std::vector<T>(n,(T)0.){}

        const T& operator ()(const int& i) const
        {
            return this->operator[](i);
        }

        T& operator ()(const int& i)
        {
            return this->operator[](i);
        }
};

#endif
