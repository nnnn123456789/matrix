#pragma once

#include <memory>
#include <functional>
#include <cstring>
#include <cmath>
#include <numeric>		
#include <tuple>

#include "locater.hpp"
#include "vector.hpp"
#include "numbers.hpp"

#pragma pack(push)
#pragma pack(4)
#pragma warning(disable:4634)




namespace yuanzm {

template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>> class matrix;

}


enum class exception{};
enum class matrix_exception{matrix_strange};

template<typename _Data, int n, typename T>
class yuanzm::matrix : public yuanzm::vector<_Data, n * n>
{
protected:
	using _MyType = matrix<_Data, n, T>;
	using _VecType = vector<_Data, n, T>;
	using _MyBase = vector<_Data, n * n>;
	constexpr static int length = n;
	constexpr static int size = length * length;
	//std::unique_ptr<_Data> head;
private:
	//_MyType& me = *this;

public:
	matrix();
	matrix(std::initializer_list<_Data> l);
	matrix(const _MyType& r);
	matrix(_MyType&& r)noexcept;
	matrix<_Data, n, T>& operator=(const matrix<_Data, n, T>& r);
	matrix<_Data, n, T>& operator=(matrix<_Data, n, T>&& r);
public:

	locater<_Data> operator[](const int n);
	const locater<const _Data> operator[](const int n) const;
	//const _Data& operator[](const int n) const;
	_Data& locate(const int x, const int y);
	const _Data& locate(const int x, const int y) const;
	matrix<_Data, n, T> operator+(const matrix<_Data, n, T>& r)const;
	matrix<_Data, n, T> operator-(const matrix<_Data, n, T>& r)const;
	///<summary>将矩阵按照所给函数初始化</summary>
	///<param name="fun">初始化矩阵使用的函数</param>
	void init(std::function<_Data(int, int)> fun);
	vector<_Data, n, T> operator*(const vector<_Data, n, T>& v) const;
	matrix<_Data, n, T> operator*(const matrix<_Data, n, T>& r) const;
	//原地转置
	void dotrans();

	///<summary>计算矩阵的转置但不改变原值</summary>
	///<return>返回转置</return>
	_MyType trans() const;

	///<seealso cref="https://www.cnblogs.com/xiaoxi666/p/6421228.html"></seealso>
	///<summary>求逆, 并返回</summary>
	_MyType LUP_solve_inverse() const;

};




template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T>& yuanzm::matrix<_Data, n, T>::operator=(const yuanzm::matrix<_Data, n, T>& r)
{
	std::memcpy(this->head.get(), r.head.get(), size * sizeof(_Data));
	return *this;
}
template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T>& yuanzm::matrix<_Data, n, T>::operator=(yuanzm::matrix<_Data, n, T>&& r)
{
	this->vector<_Data, size>::head.swap(r.vector<_Data, size>::head);
	return *this;
}


template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T> yuanzm::matrix<_Data, n, T>::operator+(const matrix<_Data, n, T>& r)const
{
	const _MyType& l = *this;
	_MyType ret;
	#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		ret[0][i] = l[0][i] + r[0][i];
	}
	return ret;
}

template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T> yuanzm::matrix<_Data, n, T>::operator-(const matrix<_Data, n, T>& r)const
{
	const _MyType& l = *this;
	_MyType ret;
	#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		ret[0][i] = l[0][i] - r[0][i];
	}
	return ret;
}


template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T> yuanzm::matrix<_Data, n, T>::operator*(const vector<_Data, n, T>& v) const
{
	vector<_Data, n, T> ret;
	#pragma omp parallel for
	for (int i = 0; i < length; i++)
	{
		ret[i] = 0;
		for (int j = 0; j < length; j++)
		{
			ret[i] += v[j] * (*this)[i][j];
		}
	}
	return ret;
}
template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T> yuanzm::matrix<_Data, n, T>::operator*(const matrix<_Data, n, T>& r) const
{
	_MyType ret;
	#pragma omp parallel for
	for (int i = 0; i < length; i++)
		for (int k = 0; k < length; k++)
		{
			ret[i][k] = 0;
			for (int j = 0; j < length; j++)
			{
				ret[i][k] += (*this)[i][j] * r[j][k];
			}
		}
	return ret;
}


template<typename _Data, int n, typename T>
yuanzm::locater<_Data> yuanzm::matrix<_Data, n, T>::operator[](const int n)
{
	auto p = this->head.get() + n * length;
	locater<_Data> ret{ p };
	return ret;
}

template<typename _Data, int n, typename T>
const yuanzm::locater<const _Data> yuanzm::matrix<_Data, n, T>::operator[](const int n) const
{
	auto p = this->head.get() + n * length;
	const locater<const _Data> ret{ p };
	return ret;
}

template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T>::matrix() : vector<_Data, size>()
{
}

template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T>::matrix(std::initializer_list<_Data> l) : matrix()
{
#pragma omp parallel for
	for (std::size_t i = 0; i < l.size() && i < size; i++)
	{
		this->head.get()[i] = l.begin()[i];
	}
}

template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T>::matrix(const _MyType& r) : vector<_Data, size>()
{
	std::memcpy(this->head.get(), r.head.get(), size * sizeof(_Data));
}

template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T>::matrix(_MyType&& r) noexcept : vector<_Data, size>(nullptr)
{
	this->head.swap(r.head);
}

template<typename _Data, int n, typename T>
_Data& yuanzm::matrix<_Data, n, T>::locate(const int x, const int y)
{
	return this->head.get()[x * length + y];
}

template<typename _Data, int n, typename T>
const _Data& yuanzm::matrix<_Data, n, T>::locate(const int x, const int y) const
{
	return this->head.get()[x * length + y];
}

///<summary>将矩阵按照所给函数初始化</summary>
///<param name="fun">初始化矩阵使用的函数</param>
template<typename _Data, int n, typename T>
void yuanzm::matrix<_Data, n, T>::init(std::function<_Data(int, int)> fun)
{
#pragma omp parallel for
	for (int i = 0; i < length; i++)
		for (int j = 0; j < length; j++)
		{
			(*this)[i][j] = fun(i, j);
		}
}

//原地转置
template<typename _Data, int n, typename T>
void yuanzm::matrix<_Data, n, T>::dotrans()
{
	yuanzm::matrix<_Data, n, T>& me = *this;
#pragma omp parallel for
	for (int i = 0; i < length; i++)
		for (int j = 0; j < i; j++)
		{
			auto temp = me[i][j];
			me[i][j] = me[j][i];
			me[j][i] = temp;
		}
}

///<summary>计算矩阵的转置但不改变原值</summary>
///<return>返回转置</return>
template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T> yuanzm::matrix<_Data, n, T>::trans() const
{
	matrix<_Data, n, T> ret;
	for (int i = 0; i < length; i++)
		for (int j = 0; j < length; j++)
			ret[i][j] = (*this)[j][i];
	return ret;
}

///<summary>求逆, 并返回</summary>
template<typename _Data, int n, typename T>
yuanzm::matrix<_Data, n, T> yuanzm::matrix<_Data, n, T>::LUP_solve_inverse() const
{
	auto LUP_Solve = [](const matrix<_Data, n, T>& L, const matrix<_Data, n, T>& U, const vector<int, n> P, const vector<_Data, n, T>& b)->vector<_Data, length>
	{
		vector<_Data, n, T> x, y;
		for (int i = 0; i < n; i++)
		{
			y[i] = b[P[i]];
			for (int j = 0; j < i; j++)
			{
				y[i] = y[i] - L[i][j] * y[j];
			}
		}
		for (int i = n - 1; i >= 0; i--)
		{
			x[i] = y[i];
			for (int j = n - 1; j > i; j--)
			{
				x[i] = x[i] - U[i][j] * x[j];
			}
			x[i] /= U[i][i];
		}
		return x;
	};
	auto LUP_Descomposition = [](matrix<_Data, n, T>& A, matrix<_Data, n, T>& L, matrix<_Data, n, T>& U, vector<int, n>& P)->void
	{
		int row = 0;
		P.init([](int i) { return i; });
		for (int i = 0; i < n - 1; i++)
		{
			long double p = 0.0;
			for (int j = i; j < n; j++)
			{
				if (std::fabsl(A[j][i]) > p)
				{
					p = std::fabsl(A[j][i]);
					row = j;
				}
			}
			if (0 == p)
			{
				throw(matrix_exception::matrix_strange);
			}
			//swap P[i] P[row]
			if (row != i)
			{
				int tmp = P[i];
				P[i] = P[row];
				P[row] = tmp;
			}
			_Data tmp2 = 0.0;
			for (int j = 0; j < n; j++)
			{
				tmp2 = A[i][j];
				A[i][j] = A[row][j];
				A[row][j] = tmp2;
			}

			long double u = A[i][i], l = 0.0;
			for (int j = i + 1; j < n; j++)
			{
				l = A[j][i] / u;
				A[j][i] = l;
				for (int k = i + 1; k < n; k++)
				{
					A[j][k] = A[j][k] - A[i][k] * l;
				}
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if (i != j)
				{
					L[i][j] = A[i][j];
				}
				else
				{
					L[i][j] = 1;
				}
			}
			for (int k = i; k < n; k++)
			{
				U[i][k] = A[i][k];
			}
		}

	};
	matrix<_Data, n, T> me_copy;
	matrix<_Data, n, T> ret; 
	vector<_Data, n, T> inv_A_each;
	vector<_Data, n, T> b;
	matrix<_Data, n, T> L;
	matrix<_Data, n, T> U;
	vector<int, n> P;
	me_copy = (*this);
	LUP_Descomposition(me_copy, L, U, P);
	for (int i = 0; i < length; i++)
	{
		b.init([]([[maybe_unused]] int x) { return 0; });
		b[i] = 1;

		inv_A_each = LUP_Solve(L, U, P, b);
		std::memcpy(ret.vector<_Data, n * n>::head.get() + i * n, &inv_A_each[0], n * sizeof(_Data));

	}
	ret.dotrans();
	return ret;
}

#pragma pack(pop)
