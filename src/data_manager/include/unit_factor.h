#pragma once
#include <map>
#include <string>
#include "types.h"

enum class UNIT
{
	LJ = 0,
	REAL,
	UNKNOWN = -1
};

std::map<std::string, UNIT unit> unit_factor_map = {
	{"LJ", UNIT::LJ},
	{"REAL", UNIT::REAL} };

template<UNIT unit>
struct UnitFactor;

template<>
struct UnitFactor<UNIT::LJ>
{
	static const rbmd::Real _kb;
	static const rbmd::Real _fmt2v;
	static const rbmd::Real _mvv2e;
	static const rbmd::Real _qqr2e;
};

const rbmd::Real UnitFactor<UNIT::LJ>::_kb = 1.0;
const rbmd::Real UnitFactor<UNIT::LJ>::_fmt2v = 1.0;
const rbmd::Real UnitFactor<UNIT::LJ>::_mvv2e = 1.0;
const rbmd::Real UnitFactor<UNIT::LJ>::_qqr2e = 1.0;

template<>
struct UnitFactor<UNIT::REAL>
{
	static const rbmd::Real _kb;
	static const rbmd::Real _fmt2v;
	static const rbmd::Real _mvv2e;
	static const rbmd::Real _qqr2e;
};

const rbmd::Real UnitFactor<UNIT::REAL>::_kb = 1.9872067 * pow(10.0, -3);
const rbmd::Real UnitFactor<UNIT::REAL>::_fmt2v = 4.186 * pow(10.0, -4);
const rbmd::Real UnitFactor<UNIT::REAL>::_mvv2e = 1.0 / (4.186 * pow(10.0, -4));
const rbmd::Real UnitFactor<UNIT::REAL>::_qqr2e = 332.06371;