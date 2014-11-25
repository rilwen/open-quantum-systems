#ifndef __TYPES_HPP
#define __TYPES_HPP

#include <boost/shared_ptr.hpp>

template <class T>
struct Typ
{
  typedef boost::shared_ptr<T> P;
};

#endif
