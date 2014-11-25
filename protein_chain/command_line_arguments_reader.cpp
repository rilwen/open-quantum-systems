#include "command_line_arguments_reader.h"

CommandLineArgumentsReader::CommandLineArgumentsReader(int argc, char* argv[])
{
	if (argc < 1) {
		throw std::domain_error("argc must be at least 1");
	}
	m_argc = argc;
	m_argv.resize(m_argc);
	for (unsigned int i = 0; i < m_argc; ++i) {
		m_argv[i] = std::string(argv[i]);
	}
	reset();
}

void CommandLineArgumentsReader::reset()
{
	m_next = m_argv.begin();
	++m_next;
}

const std::string& CommandLineArgumentsReader::first() const
{
	return m_argv.front();
}

bool CommandLineArgumentsReader::hasNext() const
{
	return m_next != m_argv.end();
}
