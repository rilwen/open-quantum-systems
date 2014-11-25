#pragma once

namespace rql
{
	namespace math
	{
		namespace sf
		{
			class NIG
			{
			public:
				NIG(double alpha, double beta, double mu, double delta);

				double alpha() const
				{
					return alpha_;
				}

				double beta() const
				{
					return beta_;
				}

				double mu() const
				{
					return mu_;
				}

				double delta() const
				{
					return delta_;
				}

				double gamma() const
				{
					return gamma_;
				}

				double pdf(double x) const;

				double mean() const
				{
					return mean_;
				}

				double variance() const
				{
					return variance_;
				}

				static double gamma(double alpha, double beta);
				static double mean(double beta, double mu, double delta, double gamma);
				static double variance(double alpha, double delta, double gamma);
			private:
				double alpha_;
				double beta_;
				double mu_;
				double delta_;
				double gamma_;
				double mean_;
				double variance_;
			};
		}
	}
}