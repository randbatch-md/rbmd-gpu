#include "unit_factor.h"

// 初始化 unit_factor_map
std::map<std::string, UNIT> unit_factor_map = { {"LJ", UNIT::LJ}, {"REAL", UNIT::REAL}, {"METAL", UNIT::METAL} };

// LJ 单位制常量定义
const rbmd::Real UnitFactor<UNIT::LJ>::_kb = 1.0;
const rbmd::Real UnitFactor<UNIT::LJ>::_fmt2v = 1.0;
const rbmd::Real UnitFactor<UNIT::LJ>::_mvv2e = 1.0;
const rbmd::Real UnitFactor<UNIT::LJ>::_qqr2e = 1.0;

// REAL 单位制常量定义
const rbmd::Real UnitFactor<UNIT::REAL>::_kb = 1.9872067 * std::pow(10.0, -3);
const rbmd::Real UnitFactor<UNIT::REAL>::_fmt2v = 4.186 * std::pow(10.0, -4);
const rbmd::Real UnitFactor<UNIT::REAL>::_mvv2e = 1.0 / (4.186 * std::pow(10.0, -4));
const rbmd::Real UnitFactor<UNIT::REAL>::_qqr2e = 332.06371;

// METAL 单位制常量定义
const rbmd::Real UnitFactor<UNIT::METAL>::_kb = 8.617343e-5;
const rbmd::Real UnitFactor<UNIT::METAL>::_fmt2v = 1.0 / 1.0364269e-4;
const rbmd::Real UnitFactor<UNIT::METAL>::_mvv2e = 1.0364269e-4;
const rbmd::Real UnitFactor<UNIT::METAL>::_qqr2e = 14.399645;
