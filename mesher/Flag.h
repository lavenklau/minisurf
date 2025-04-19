#pragma once
#include "Config.h"

BEGIN_MINSURF_NAMESPACE

struct FlagBase {
	int _flag = 0;
	enum flag_t {
		xminBoundary = 1, xmaxBoundary = 0b10,
		yminBoundary = 0b100, ymaxBoundary = 0b1000,
		zminBoundary = 0b10000, zmaxBoundary = 0b100000,
		periodMask = 0b111111,
		minPeriodMask = 0b010101,
		maxPeriodMask = 0b101010,
		D_BOUNDARY = 0b1000000,
		N_BOUNDARY = 0b10000000,
		ChartMask =  0b11100000000,
	};
	bool is_period_boundary(void) {
		return _flag & flag_t::periodMask;
	}
	bool is_period_boundary(int axis) {
		return _flag & (0b11 << (axis * 2));
	}
	bool is_period_edge(void) const { return (bool(_flag & 0b11) + bool(_flag & 0b1100) + bool(_flag & 0b110000)) >= 2; }
	int period_type_id(void) const { return (_flag & 0b11) + ((_flag & 0b1100) >> 2) * 3 + ((_flag & 0b110000) >> 4) * 9; }
	void set_period_id_type(int id) { set_period_boundary(id % 3, id / 3 % 3, id / 9); }
	bool is_min_period(void) const { return _flag & flag_t::minPeriodMask; }
	bool is_min_period(int axis) const { return _flag & (1 << (axis * 2)); }
	bool is_max_period(void) const { return _flag & flag_t::maxPeriodMask; }
	bool is_max_period(int axis) const { return _flag & (2 << (axis * 2)); }
	bool is_d_boundary(void) const { return _flag & D_BOUNDARY; }
	bool is_n_boundary(void) const { return _flag & N_BOUNDARY; }
	void set_d_boundary(void) { _flag = _flag | D_BOUNDARY; }
	void set_n_boundary(void) { _flag = _flag | N_BOUNDARY; }
	int getChart(void) const { return (_flag & flag_t::ChartMask) >> 8; }
	void setChart(int x, int y, int z) {
		_flag &= ~flag_t::ChartMask;
		_flag |= (bool(x) << 8) | (bool(y) << 9) | (bool(z) << 10);
	}
	template<typename Scalar>
	void set_period_boundary(Scalar border, Scalar px, Scalar py, Scalar pz) {
		_flag = (_flag & ~periodMask) |
			((px > border) << 1) | ((px < -border) << 0) |
			((py > border) << 3) | ((py < -border) << 2) |
			((pz > border) << 5) | ((pz < -border) << 4);
	}
	void set_period_boundary(int x, int y, int z) {
		_flag &= ~periodMask;
		_flag |= (x | (y << 2) | (z << 4)) & periodMask;
	}
	void set_min_boundary(int axis) { _flag |= 1 << (axis * 2); }
	void set_max_boundary(int axis) { _flag |= 2 << (axis * 2); }
	void set_period_boundary(int xyz) {
		_flag &= ~periodMask;
		_flag |= xyz & periodMask;
	}
	int getPeriodFlag(void) {
		return _flag & flag_t::periodMask;
	}
};

struct VertexFlag : FlagBase { };
struct EdgeFlag : FlagBase {};

END_MINSURF_NAMESPACE
