#include "tictoc.h"

std::map<std::string, double> tictoc::Record::_table;

std::stack<std::unique_ptr<tictoc::live>> tictoc::_timerstack;

bool tictoc::_mode_accumulation = false;

tictoc::void_buf tictoc::_voidBuf;

std::ostream  tictoc::_silent_stream(&_voidBuf);

std::chrono::steady_clock::time_point tictoc::getTag(void)
{
	return std::chrono::steady_clock::now();
}

double tictoc::get_record(const std::string& rec_name)
{
	auto it = _record._table.find(rec_name);
	if (it != _record._table.end()) {
		return it->second;
	}
	else {
		return 0;
	}
}

std::map<std::string, double> tictoc::clear_record(void)
{
	std::map<std::string, double> oldmap = _record._table;
	_record._table.clear();
	return oldmap;
}
