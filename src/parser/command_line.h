#pragma once

class CommandLine
{
public:
	CommandLine(int argc, char* argv[]);
	virtual ~CommandLine() = default;

public:
	bool IsRunApplication() { return _is_run_application; }

private:
	bool _is_run_application;
};