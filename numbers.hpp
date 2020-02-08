#pragma once
//#include <concepts>
#include <type_traits>

namespace yuanzm
{
template<typename T>
using enable_if_floating_point = ::std::enable_if<::std::is_floating_point_v<double>>;
template<typename T>
using enable_if_floating_point_t = typename enable_if_floating_point<T>::type;
}
namespace yuanzm::numbers {

template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T e_v =          2.7182818284590452353602874713527;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T log2e_v =      1.4426950408889634073599246810019;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T log10e_v =     .43429448190325182765112891891661;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T pi_v =         3.1415926535897932384626433832795;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T inv_pi_v =     .31830988618379067153776752674503;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T inv_sqrtpi_v = .56418958354775628694807945156077;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T ln2_v =        .69314718055994530941723212145818;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T ln10_v =       2.3025850929940456840179914546844;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T sqrt2_v =      1.4142135623730950488016887242097;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T sqrt3_v =      1.7320508075688772935274463415059;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T inv_sqrt2_v =  .70710678118654752440084436210485;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T inv_sqrt3_v =  .57735026918962576450914878050196;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T egamma_v =     .57721566490153286060651209008240;
template<class T, typename = enable_if_floating_point_t<T>> 
constexpr T phi_v =        1.6180339887498948482045868343656;

constexpr double e = e_v<double>;
constexpr double log2e = log2e_v<double>;
constexpr double log10e = log10e_v<double>;
constexpr double pi = pi_v<double>;
constexpr double inv_pi = inv_pi_v<double>;
constexpr double inv_sqrtpi = inv_sqrtpi_v<double>;
constexpr double ln2 = ln2_v<double>;
constexpr double ln10 = ln10_v<double>;
constexpr double sqrt2 = sqrt2_v<double>;
constexpr double sqrt3 = sqrt3_v<double>;
constexpr double inv_sqrt2 = inv_sqrt2_v<double>;
constexpr double inv_sqrt3 = inv_sqrt3_v<double>;
constexpr double egamma = egamma_v<double>;
constexpr double phi = phi_v<double>;

}

