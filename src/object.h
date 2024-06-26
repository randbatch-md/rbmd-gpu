#pragma once
#include <memory>
#include <iostream>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

class Object
{
public:
	Object()
	{
		_console = spdlog::get("console");
		if (!_console)
		{
			_console = spdlog::stdout_color_mt("console");
		}
	}

	virtual ~Object() = default;

protected:
	std::shared_ptr<spdlog::logger> _console;
};