#ifndef __RQL_MATH_PIMPL_H
#define __RQL_MATH_PIMPL_H

#include <memory>
#include <stdexcept>

namespace rql {
	template <class Impl>
	class Pimpl
	{
	public:
		Pimpl();
		Pimpl(std::shared_ptr<Impl> impl);
		Pimpl(const Pimpl<Impl>& other);
		Pimpl<Impl>& operator=(const Pimpl<Impl>& that);		
		void swap(Pimpl<Impl>& other);
		bool operator==(const Pimpl<Impl>& other) const;
		bool operator!=(const Pimpl<Impl>& other) const;
	protected:
		Impl& impl();
		const Impl& impl() const;
	private:
		std::shared_ptr<Impl> m_impl;
	};

	template <class Impl> Pimpl<Impl>::Pimpl()
		: m_impl(std::shared_ptr<Impl>())
	{
	}

	template <class Impl> Pimpl<Impl>::Pimpl(std::shared_ptr<Impl> impl)
		: m_impl(impl)
	{
		if (!m_impl)
			throw std::runtime_error("Pimpl: null implementation");
	}

	template <class Impl> Pimpl<Impl>::Pimpl(const Pimpl<Impl>& other)
		: m_impl(other.m_impl)
	{
	}

	template <class Impl> Pimpl<Impl>& Pimpl<Impl>::operator=(const Pimpl<Impl>& that)
	{
		m_impl = that.m_impl;
		return *this;
	}

	template <class Impl> void Pimpl<Impl>::swap(Pimpl<Impl>& other)
	{
		m_impl.swap(other.m_impl);
	}

	template <class Impl> Impl& Pimpl<Impl>::impl()
	{
		return *m_impl;
	}

	template <class Impl> const Impl& Pimpl<Impl>::impl() const
	{
		return *m_impl;
	}

	template <class Impl> bool Pimpl<Impl>::operator==(const Pimpl<Impl>& other) const
	{
		return m_impl == other.m_impl;
	}

	template <class Impl> bool Pimpl<Impl>::operator!=(const Pimpl<Impl>& other) const
	{
		return m_impl != other.m_impl;
	}
}

#endif // __RQL_MATH_PIMPL_H