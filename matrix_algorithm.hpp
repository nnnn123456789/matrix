#pragma once
#include "vector.hpp"
#include "matrix.hpp"
#define MY_DEBUG
#ifdef MY_DEBUG
#include <iostream>
#endif
namespace yuanzm {
constexpr int iter_stop = 10000;
constexpr long double precision = 1e-6;

template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>>
[[nodiscard]] std::unique_ptr<yuanzm::vector<_Data, n>> solve_catch_ptr(const yuanzm::vector<_Data, n - 1>& _d, const yuanzm::vector<_Data, n>& _c, const yuanzm::vector<_Data, n - 1>& _u, const yuanzm::vector<_Data, n >& _b);

template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>>
[[nodiscard]] yuanzm::vector<_Data, n> solve_catch(const yuanzm::vector<_Data, n - 1>& _d, const yuanzm::vector<_Data, n>& _c, const yuanzm::vector<_Data, n - 1>& _u, const yuanzm::vector<_Data, n >& _b);

template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>>
[[nodiscard]] std::tuple<vector<_Data, n>, int> Jacobi_method(const matrix<_Data, n>& A, const vector<_Data, n>& b);

template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>>
[[nodiscard]] std::tuple<vector<_Data, n>, int> Gauss_Seidel_method(const matrix<_Data, n>& A, const vector<_Data, n>& b);

template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>>
[[nodiscard]] std::tuple<vector<_Data, n>, int> SOR_method(const matrix<_Data, n>& A, const vector<_Data, n>& b, const double omega = 1);

}

template<typename _Data, int n, typename>
std::unique_ptr<yuanzm::vector<_Data, n>> yuanzm::solve_catch_ptr(
	const yuanzm::vector<_Data, n - 1>& _d, 
	const yuanzm::vector<_Data, n>& _c, 
	const yuanzm::vector<_Data, n - 1>& _u, 
	const yuanzm::vector<_Data, n >& _b)
{
	std::unique_ptr<yuanzm::vector<_Data, n - 1>> pd(new yuanzm::vector<_Data, n - 1>);
	std::unique_ptr<yuanzm::vector<_Data, n>> pc(new yuanzm::vector<_Data, n>);
	std::unique_ptr<yuanzm::vector<_Data, n - 1>> pu(new yuanzm::vector<_Data, n - 1>);
	std::unique_ptr<yuanzm::vector<_Data, n>> pb(new yuanzm::vector<_Data, n>);
	auto& d = *pd; d = _d;
	auto& c = *pc; c = _c;
	auto& u = *pu; u = _u;
	auto& b = *pb; b = _b;
	std::unique_ptr<yuanzm::vector<_Data, n>> px(new yuanzm::vector<_Data, n>);
	auto& x = *px;
	for (int i = 1; i < n; i++)
	{
		d[i - 1] /= c[i - 1];
		c[i] -= u[i - 1] * d[i - 1];
		b[i] -= b[i - 1] * d[i - 1];
	}
	x[n - 1] = b[n - 1] / c[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (b[i] - u[i] * x[i + 1]) / c[i];
	}
	return std::move(px);
}

///<return> 1: answer, 2: iter-count </return>
template<typename _Data, int n, typename T>
std::tuple<yuanzm::vector<_Data, n>, int> yuanzm::Jacobi_method(const yuanzm::matrix<_Data, n>& A, const yuanzm::vector<_Data, n>& b)
{
	vector<_Data, n> last, current;

	last = b;
	int t = 0;
	while (++t < iter_stop)
	{
		for (int i = 0; i < n; i++)
		{
			_Data s = 0;
			for (int j = 0; j < n; j++)
			{
				if (j == i)
					continue;
				else
					s += A[i][j] * last[j];
			}
			current[i] = (b[i] - s) / A[i][i];
		}
		if (vector<_Data, n>::supnorm(last - current) < precision)
			break;
		else
			last = current;
	}
	return std::make_tuple(current, t);
}


///<return> 1: answer, 2: iter-count </return>
template<typename _Data, int n, typename T>
std::tuple<yuanzm::vector<_Data, n>, int> yuanzm::Gauss_Seidel_method(const yuanzm::matrix<_Data, n>& A, const yuanzm::vector<_Data, n>& b)
{
	vector<_Data, n> last, current;

	last = b;
	int t = 0;
	while (t++ < iter_stop)
	{
		for (int i = 0; i < n; i++)
		{
			_Data s = 0;
			for (int j = 0; j < n; j++)
			{
				if (j == i)
					continue;
				else
					s += A[i][j] * current[j];
			}
			current[i] = (b[i] - s) / A[i][i];
		}
		auto temp = vector<_Data, n>::supnorm(last - current);
#ifdef MY_DEBUG
		//std::cerr << current << std::endl;
		std::cerr << temp << std::endl;
#endif
		if (temp < precision)
			break;
		else
			last = current;
	}
	return std::make_tuple(current, t);
}

///<return> 1: answer, 2: iter-count </return>
template<typename _Data, int n, typename T>
std::tuple<yuanzm::vector<_Data, n>, int> yuanzm::SOR_method(const yuanzm::matrix<_Data, n>& A, const yuanzm::vector<_Data, n>& b, const double omega)
{
	vector<_Data, n> last, current;

	last = b;
	int t = 0;
	while (t++ < iter_stop)
	{
		for (int i = 0; i < n; i++)
		{
			_Data s = 0;
			//#pragma omp parallel for
			for (int j = 0; j < n; j++)
			{
				if (j == i)
					continue;
				else
					s += A[i][j] * current[j];
			}
			current[i] = (b[i] - s) / A[i][i];
			current[i] = (1 - omega) * last[i] + omega * current[i];
		}
		auto temp = vector<_Data, n>::supnorm(last - current);
#ifdef MY_DEBUG
		//std::cerr << current << std::endl;
		std::cerr << temp << std::endl;
#endif
		if (temp < precision)
			break;
		else
			last = current;
	}
	return std::make_tuple(current, t);
} 


template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n> yuanzm::solve_catch(const yuanzm::vector<_Data, n - 1>& _d, const yuanzm::vector<_Data, n>& _c, const yuanzm::vector<_Data, n - 1>& _u, const yuanzm::vector<_Data, n>& _b)
{
	auto d = _d;
	auto c = _c;
	auto&u = _u;
	auto b = _b;
	decltype(c) x;
	for (int i = 1; i < n; i++)
	{
		d[i - 1] /= c[i - 1];
		c[i] -= u[i - 1] * d[i - 1];
		b[i] -= b[i - 1] * d[i - 1];
	}
	x[n - 1] = b[n - 1] / c[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (b[i] - u[i] * x[i + 1]) / c[i];
	}
	return std::move(x);
}
