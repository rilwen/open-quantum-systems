#ifndef __COMMAND_LINE_ARGUMENTS_READER_H
#define __COMMAND_LINE_ARGUMENTS_READER_H

#include <vector>
#include <string>
#include <stdexcept>
#include "utils.h"
#include "core.h"

class CommandLineArgumentsReader
{
public:
	PROTEIN_CHAIN_API CommandLineArgumentsReader(int argc, char* argv[]);
	void reset();
	const std::string& first() const;
	PROTEIN_CHAIN_API bool hasNext() const;
	template <class T> T pop();
	template <class T> void pop(T& dest);
private:
	unsigned int m_argc;
	std::vector<std::string> m_argv;
	std::vector<std::string>::const_iterator m_next;
};

template <class T> T CommandLineArgumentsReader::pop()
{
	T dest;
	pop(dest);
	return dest;
}

template <class T> void CommandLineArgumentsReader::pop(T& dest)
{
	if (hasNext()) {
		dest = strConvert<T>(*m_next);
		++m_next;
	} else {
		throw std::runtime_error("Run out of arguments");
	}
}


#endif // __COMMAND_LINE_ARGUMENTS_READER_H
