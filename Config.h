#pragma once

#include <stdio.h>
#include "numeric"
#include "limits"

#define BEGIN_MINSURF_NAMESPACE namespace msf {

#define END_MINSURF_NAMESPACE }

BEGIN_MINSURF_NAMESPACE

using Real = double;
constexpr Real min_real = -std::numeric_limits<Real>::max();
constexpr Real max_real = std::numeric_limits<Real>::max();
constexpr Real min_abs = std::numeric_limits<Real>::min();
constexpr Real eps_real = 1e-7;

END_MINSURF_NAMESPACE


#ifdef __linux__
template<int N, typename... Args>
void sprintf_s(char(&_Buffer)[N], Args... args) { snprintf(_Buffer, N, args...); }
#elif defined(_WIN32)
#else
#error Only Windows and Linux are supported!
#endif


#if _WIN32
#define SELECTANY __declspec(selectany) 
#elif __linux__
#define SELECTANY __attribute__((weak))
#endif
