#pragma once
#include "numbers.hpp"

namespace yuanzm
{



template<typename T1, typename T2>
auto constexpr max(T1 l, T2 r) ->decltype(T1{} +T2{})
{
	return l > r ? l : r;
}

template<typename T1, typename ... types>
auto constexpr max(T1 l, types... paras) ->decltype(T1{} +max(paras...))
{
	decltype(auto) t = max(paras...);
	return l > t ? l : t;
}

//template<typename T>
////requires {T::begin(); T::end(); }
//auto constexpr max(T iter) ->decltype(*iter.begin())
//{
//	auto ret = *(iter.begin());
//	for (auto i : iter)
//	{
//		ret = max(ret, i);
//	}
//	return ret;
//}


template<typename T>
T abs(T x)
{
	return x > T{} ? x : -x;
}



}
