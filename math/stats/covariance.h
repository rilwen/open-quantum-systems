#ifndef __RQL_MATHS_COVARIANCE_H
#define __RQL_MATHS_COVARIANCE_H

#include "../average.h"
#include "../MathCore.h"
#include <cassert>
#include <vector>
#include <boost/multi_array.hpp>
#include <Eigen/Core>

namespace rql {
	namespace math {
		namespace stats {

			class Covariance
			{
			public:
				typedef Eigen::MatrixXd matrix_type;
				RQL_MATH_API Covariance(unsigned int dim);
				RQL_MATH_API unsigned int dim() const { return m_dim; }
				template <class Iter1, class Iter2> void update(const Iter1 begin, const Iter2 end);
				RQL_MATH_API matrix_type covariance() const;
			private:
				Covariance& operator=(const Covariance&); // not implemented
				const unsigned int m_dim;
				std::vector<Average<double> > m_avg1;
				boost::multi_array<Average<double>, 2> m_avg2;
			};

			template <class Iter1, class Iter2> void Covariance::update(const Iter1 begin, const Iter2 end)
			{
				std::vector<Average<double> >::iterator a1it = m_avg1.begin();
				unsigned int r = 0;
				for (Iter1 i = begin; i != end; ++i) {
					const double xi = *i;
					assert( a1it != m_avg1.end() );
					a1it->update(xi);
					assert( r < m_dim );
					m_avg2[r][r].update(xi*xi);
					unsigned int c = r + 1;
					for (Iter1 j = i+1; j != end; ++j) {
						assert( c < m_dim );
						m_avg2[r][c].update(xi*(*j));
						++c;
					}
					++r;
					++a1it;
				}
			}

		}
	}
}
	

#endif // __RQL_MATHS_COVARIANCE_H
