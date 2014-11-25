#include "chaos_calculator.h"
#include <stdexcept>

ChaosCalculator::ChaosCalculator(unsigned int size)
	: m_fft(size, FFTW_FORWARD, FFTW_MEASURE), m_deltas(size), m_spectrum(size), m_averages(size)
{
	m_deltas[0] = 0;	
}


