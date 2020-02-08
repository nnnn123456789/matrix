#pragma once
namespace yuanzm
{
template<typename _Data> class locater;
}

template<typename _Data>
class yuanzm::locater
{
private:
	_Data* s;
public:
	locater(_Data* p);
	_Data& operator[](const int n);
	const _Data& operator[](const int n) const;
};

template<typename _Data>
yuanzm::locater<_Data>::locater(_Data* p) : s(p) {}

template<typename _Data>
_Data& yuanzm::locater<_Data>::operator[](const int n)
{
	return s[n];
}

template<typename _Data>
const _Data& yuanzm::locater<_Data>::operator[](const int n) const
{
	return s[n];
}
