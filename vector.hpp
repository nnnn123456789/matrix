#pragma once

#include "numbers.hpp"
#include <algorithm>
namespace yuanzm
{
template<typename _Data, int n, typename = enable_if_floating_point_t<_Data>> class vector;
}


template<typename _Data, int n, typename T>
class yuanzm::vector
{

protected:
	using _MyType = vector<_Data, n, T>;
	constexpr static int length = n;
	std::unique_ptr<_Data> head;
private:
	//_MyType& me = *this;
public:
	vector();
	vector(std::initializer_list<_Data> l);
	vector(const _MyType& r);
	vector(_MyType&& r)noexcept;
	vector(std::nullptr_t);
	vector<_Data, n, T>& operator=(const vector<_Data, n, T>& r);
	vector<_Data, n, T>& operator=(vector<_Data, n, T>&& r)noexcept;
	_Data& operator[](const int n);
	_Data& locate(const int n);
	const _Data& operator[](const int n) const;
	const _Data& locate(const int n) const;
	vector<_Data, n, T> operator+(const vector<_Data, n, T>& r)const;
	vector<_Data, n, T> operator-(const vector<_Data, n, T>& r)const;
	void init(std::function<_Data(int)> fun);
	_Data* begin();
	_Data* end();
	const _Data* begin()const;
	const _Data* end()const;
	_Data supnorm() const;
	_Data norm(std::function<_Data(_Data, _Data)> add) const;
	static _Data supnorm(const vector<_Data, n, T>& vec);
	static _Data norm(const vector<_Data, n, T>& vec, std::function<_Data(_Data, _Data)> add);
};


template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>::vector() : head(new _Data[length]) {
	//std::for_each(head.get(), head.get() + length, [](_Data& x) {x = _Data{}; });
	std::for_each_n(head.get(), length, [](_Data& x) {x = _Data{}; });

}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>::vector(std::initializer_list<_Data> l) : head(new _Data[length])
{
#pragma omp parallel for
	for (std::size_t i = 0; i < l.size() && i < length; i++)
	{
		head.get()[i] = l.begin()[i];
	}
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>::vector(const _MyType& r) : head(new _Data[length])
{
	std::memcpy(this->head.get(), r.head.get(), length * sizeof(_Data));
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>::vector(_MyType&& r) noexcept : head(nullptr)
{
	this->head.swap(r.head);
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>::vector(std::nullptr_t) :head(nullptr)
{

}



template<typename _Data, int n, typename T>
_Data& yuanzm::vector<_Data, n, T>::locate(const int n)
{
	return head.get()[n];
}

template<typename _Data, int n, typename T>
const _Data& yuanzm::vector<_Data, n, T>::locate(const int n) const
{
	return head.get()[n];
}

template<typename _Data, int n, typename T>
void yuanzm::vector<_Data, n, T>::init(std::function<_Data(int)> fun)
{
#pragma omp parallel for
	for (int i = 0; i < length; i++)
	{
		(*this)[i] = fun(i);
	}
}

template<typename _Data, int n, typename T>
_Data* yuanzm::vector<_Data, n, T>::begin()
{
	return head.get();
}

template<typename _Data, int n, typename T>
_Data* yuanzm::vector<_Data, n, T>::end()
{
	return head.get() + length;
}

template<typename _Data, int n, typename T>
const _Data* yuanzm::vector<_Data, n, T>::begin()const
{
	return head.get();
}

template<typename _Data, int n, typename T>
const _Data* yuanzm::vector<_Data, n, T>::end()const
{
	return head.get() + length;
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>& yuanzm::vector<_Data, n, T>::operator=(vector<_Data, n, T>&& r)noexcept
{
	this->head.swap(r.head);
	return *this;
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T>& yuanzm::vector<_Data, n, T>::operator=(const yuanzm::vector<_Data, n, T>& r)
{
	std::memcpy(this->head.get(), r.head.get(), length * sizeof(_Data));
	return *this;
}

template<typename _Data, int n, typename T>
const _Data& yuanzm::vector<_Data, n, T>::operator[](const int n) const
{
	return head.get()[n];
}

template<typename _Data, int n, typename T>
_Data& yuanzm::vector<_Data, n, T>::operator[](const int n)
{
	return head.get()[n];
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T> yuanzm::vector<_Data, n, T>::operator+(const vector<_Data, n, T>& r)const
{
	const _MyType& l = *this;
	_MyType ret;
#pragma omp parallel for
	for (int i = 0; i < length; i++)
	{
		ret[i] = l[i] + r[i];
	}
	return ret;
}

template<typename _Data, int n, typename T>
yuanzm::vector<_Data, n, T> yuanzm::vector<_Data, n, T>::operator-(const vector<_Data, n, T>& r) const
{
	const _MyType& l = *this;
	_MyType ret;
#pragma omp parallel for
	for (int i = 0; i < length; i++)
	{
		ret[i] = l[i] - r[i];
	}
	return ret;
}



template<typename _Data, int n, typename T>
_Data yuanzm::vector<_Data, n, T>::norm(std::function<_Data(_Data, _Data)> add) const
{
	_Data s = std::accumulate(this->begin(), this->end(), (_Data)0, add);
	return s;
}

template<typename _Data, int n, typename T>
_Data yuanzm::vector<_Data, n, T>::supnorm() const
{
	auto f = [](_Data l, _Data r)->_Data
	{
		return yuanzm::max(yuanzm::abs(l), yuanzm::abs(r));
	};
	return this->norm(f);
	//return yuanzm::max(*this);
}

template<typename _Data, int n, typename T>
_Data yuanzm::vector<_Data, n, T>::supnorm(const vector<_Data, n, T>& vec)
{
	return vec.supnorm();
}

template<typename _Data, int n, typename T>
_Data yuanzm::vector<_Data, n, T>::norm(const vector<_Data, n, T>& vec, std::function<_Data(_Data, _Data)> add)
{
	return vec.norm(add);
}

template<typename _char, typename _trails, typename _Data, int n>
std::basic_ostream<_char, _trails>& operator<<(std::basic_ostream<_char, _trails>& os,const yuanzm::vector<_Data, n>& v)
{
	for (auto i : v)
	{
		os << i << '\t';
	}
	return os;
}
