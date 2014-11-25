#pragma once
#pragma warning(disable:4996)

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <ostream>
#include <boost/noncopyable.hpp>
#include <boost/scoped_array.hpp>
#include "MathUtils.h"

namespace rql { namespace math {

	template <class T> class Jagged2DArray;

	template <class T>
	class Jagged2DArrayRowRef: public boost::noncopyable
	{
	public:
		Jagged2DArrayRowRef(const Jagged2DArrayRowRef<T>& other);

		// Empty reference.
		Jagged2DArrayRowRef();

		~Jagged2DArrayRowRef();

		size_t size() const;
		T& operator[](const size_t i);
		const T& operator[](const size_t i) const;

		typedef T* iterator;
		typedef const T* const_iterator;
		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;

		/** Copy the reference data from other object */
		void reset(const Jagged2DArrayRowRef<T>& other);
		friend class Jagged2DArray<T>;
	private:
		Jagged2DArrayRowRef(T* const data, const size_t size);
		T* data_;
		size_t size_;
	};

	template <class T>
	Jagged2DArrayRowRef<T>::Jagged2DArrayRowRef(T* const data, const size_t size)
	: data_(data), size_(size)
	{
	}

	template <class T>
	Jagged2DArrayRowRef<T>::Jagged2DArrayRowRef()
	: data_(0), size_(0)
	{
	}

	template <class T>
	Jagged2DArrayRowRef<T>::Jagged2DArrayRowRef(const Jagged2DArrayRowRef<T>& other)
	: data_(other.data_), size_(other.size_)
	{
	}

	template <class T>
	Jagged2DArrayRowRef<T>::~Jagged2DArrayRowRef()
	{
	}

	template <class T>
	inline size_t Jagged2DArrayRowRef<T>::size() const
	{
		return size_;
	}

	template <class T>
	inline T& Jagged2DArrayRowRef<T>::operator[](const size_t i)
	{
	#ifndef NDEBUG
		if (i >= size_) {
			throw std::out_of_range("Row element index out of range");
		}
	#endif
		return data_[i];
	}

	template <class T>
	inline const T& Jagged2DArrayRowRef<T>::operator[](const size_t i) const
	{
	#ifndef NDEBUG
		if (i >= size_) {
			throw std::out_of_range("Row element index out of range");
		}
	#endif
		return data_[i];
	}

	template <class T>
	inline T* Jagged2DArrayRowRef<T>::begin()
	{
		return data_;
	}

	template <class T>
	inline T* Jagged2DArrayRowRef<T>::end()
	{
		return data_ + size_;
	}

	template <class T>
	inline const T* Jagged2DArrayRowRef<T>::begin() const
	{
		return data_;
	}

	template <class T>
	inline const T* Jagged2DArrayRowRef<T>::end() const
	{
		return data_ + size_;
	}

	template <class T>
	void Jagged2DArrayRowRef<T>::reset(const Jagged2DArrayRowRef<T>& other)
	{
		data_ = other.data_;
		size_ = other.size_;
	}	
	
	template <class T>
	class Jagged2DArray
	{
	public:
		Jagged2DArray(const size_t nbr_rows, const size_t nbr_cols, const T& elem );
		Jagged2DArray(const size_t nbr_rows, const size_t nbr_cols);
		//! Create a lower triangular matrix
		//! @param N number rows
		//! @param include_diagonal Does the matrix include diagonal entries (i == j) or just those below the diagonal (i > j)?
		Jagged2DArray(const size_t N, const bool include_diagonal);
		Jagged2DArray(const size_t N, const bool include_diagonal, const T& elem);
		Jagged2DArray();
		Jagged2DArray(const Jagged2DArray<T>& other);
		////! Copy dimensions from other jagged array
		//template <class S> Jagged2DArray(const Jagged2DArray<S>& other);

		/**
		 * Copy the row sizes from the iterator.
		 * @tparam Iter Const Iterator.
		 */
		template <class Iter>
		Jagged2DArray(Iter begin, Iter end);

		/**
		Copy the data from the old container.
		@tparam C 2D container with begin(), end() and size() on each level.
		template <class C>
		*/
		template <class C>
		Jagged2DArray(const C& data);

		/** Does not throw */
		void swap(Jagged2DArray<T>& other);

		Jagged2DArray<T>& operator=(const Jagged2DArray<T>& other);				

		~Jagged2DArray();

		/** Number of rows */
		size_t size() const;

		/** Checks indices in debug mode
		 * @throws std::out_of_range if index out of range
		*/
		size_t row_size(const size_t row) const;

		/** Total number of elements */
		size_t nbr_elements() const;

		/** Checks indices in debug mode
		 * @throws std::out_of_range if index out of range
		*/
		T& operator()(const size_t row, const size_t row_elem);
		/** Checks indices in debug mode
		 * @throws std::out_of_range if index out of range
		*/
		const T& operator()(const size_t row, const size_t row_elem) const;

		//! Fill with constant element
		void fill(const T& elem);

		Jagged2DArrayRowRef<T> operator[](const size_t row_idx);
		const Jagged2DArrayRowRef<T> operator[](const size_t row_idx) const;

		Jagged2DArrayRowRef<T> flat_form();
		const Jagged2DArrayRowRef<T> flat_form() const;
		typedef typename Jagged2DArrayRowRef<T>::iterator flat_iterator;
		flat_iterator flat_begin();
		flat_iterator flat_end();
		typedef typename Jagged2DArrayRowRef<T>::const_iterator const_flat_iterator;
		const_flat_iterator flat_begin() const;
		const_flat_iterator flat_end() const;

		// friends
		friend Jagged2DArray<T>& operator+=(Jagged2DArray<T>& left, const Jagged2DArray<T>& right) {
		  assert(left.size_ == right.size_);
		  const T* right_ptr = right.data_.get();
		  T* const left_end = left.data_.get() + left.size_;
		  for (T* it = left.data_.get(); it != left_end; ++it) {
			  assert( right_ptr != right.data_.get() + right.size_ );
			  *it += *right_ptr;
			  ++right_ptr;
		  }
		  return left;
		}
		
		friend Jagged2DArray<T>& operator-=(Jagged2DArray<T>& left, const Jagged2DArray<T>& right) {
		  assert(left.size_ == right.size_);
		  const T* right_ptr = right.data_.get();
		  T* const left_end = left.data_.get() + left.size_;
		  for (T* it = left.data_.get(); it != left_end; ++it) {
			  assert( right_ptr != right.data_.get() + right.size_ );
			  *it -= *right_ptr;
			  ++right_ptr;
		  }
		  return left;
		}
		
		friend Jagged2DArray<T>& operator*=(Jagged2DArray<T>& left, const T& right)
		{
		  T* const left_end = left.data_.get() + left.size_;
		  for (T* it = left.data_.get(); it != left_end; ++it) {
			  *it *= right;
		  }
		  return left;
		}
		
		friend Jagged2DArray<T>& operator/=(Jagged2DArray<T>& left, const T& right)
		{
		  T* const left_end = left.data_.get() + left.size_;
		  for (T* it = left.data_.get(); it != left_end; ++it) {
			  *it /= right;
		  }
		  return left;
		}
	private:
		void init_rectangular(const size_t nbr_rows, const size_t nbr_cols);
		void init_triangular(const size_t nbr_rows, const bool include_diagonal);
		static size_t triangular_size(size_t N, bool include_diagonal);
		/** Number of rows */
		size_t size_;
		boost::scoped_array<size_t> row_starts_;
		boost::scoped_array<T> data_;
	};

	template <class T>
	Jagged2DArray<T>::Jagged2DArray(const size_t nbr_rows, const size_t nbr_cols, const T& elem)
	{
		init_rectangular(nbr_rows, nbr_cols);
		fill(elem);
	}

	template <class T>
	Jagged2DArray<T>::Jagged2DArray(const size_t nbr_rows, const size_t nbr_cols)
	{
		init_rectangular(nbr_rows, nbr_cols);
	}

	template <class T>
	Jagged2DArray<T>::Jagged2DArray(const size_t N, const bool include_diagonal)
	{
		init_triangular(N, include_diagonal);
	}

	template <class T>
	Jagged2DArray<T>::Jagged2DArray(const size_t N, const bool include_diagonal, const T& elem)
	{
		init_triangular(N, include_diagonal);
		fill(elem);
	}

	template <class T>
	Jagged2DArray<T>::Jagged2DArray()
	{
		init_rectangular(0u, 0u);
	}

	template <class T>
	Jagged2DArray<T>::Jagged2DArray(const Jagged2DArray<T>& other)
	: size_(other.size_), row_starts_(new size_t[other.size_ + 1]), data_(new T[other.row_starts_[other.size_]])
	{
		std::copy(other.row_starts_.get(), other.row_starts_.get() + size_ + 1, row_starts_.get());
		std::copy(other.data_.get(), other.data_.get() + row_starts_[size_], data_.get());
	}

	/*template <class T>
	template <class S> Jagged2DArray<T>::Jagged2DArray(const Jagged2DArray<S>& other)
		: size_(other.size_), row_starts_(new size_t[other.size_ + 1]), data_(new T[other.row_starts_[other.size_]])
	{
		std::copy(other.row_starts_.get(), other.row_starts_.get() + size_ + 1, row_starts_.get());
	}*/

	template <class T>
	template <class Iter>
	Jagged2DArray<T>::Jagged2DArray(Iter begin, Iter end)
	: size_(0), row_starts_(0), data_(0)
	{
		// count rows
		size_t idx = 0;
		for (Iter i = begin; i != end; ++i) {
			++idx;
		}
		size_ = idx;
		if (size_ > 0) {
			row_starts_.reset(new size_t[size_ + 1]);
			size_t total_data_size = 0;
			idx = 0;
			for (Iter i = begin; i != end; ++i) {
				size_t row_size = *i;
				row_starts_[idx] = total_data_size;
				total_data_size += row_size;
				++idx;
			}
			row_starts_[idx] = total_data_size;
			if (total_data_size > 0) {
				data_.reset(new T[total_data_size]);
			}
		}
	}

	template <class T>
	template <class C>
	Jagged2DArray<T>::Jagged2DArray(const C& data)
	: size_(data.size()), row_starts_(new size_t[size_ + 1]), data_(0)
	{
		size_t total_data_size = 0;
		size_t idx = 0;
		for (typename C::const_iterator i = data.begin(); i != data.end(); ++i) {
			const size_t row_size = i->size();
			row_starts_[idx] = total_data_size;
			total_data_size += row_size;
			++idx;
		}
		row_starts_[idx] = total_data_size;
		if (total_data_size > 0) {
			data_.reset(new T[total_data_size]);
			idx = 0;
			for (typename C::const_iterator i = data.begin(); i != data.end(); ++i) {
				std::copy(i->begin(), i->end(), data_.get() + row_starts_[idx]);
				++idx;
			}
		}
	}

	template <class T>
	Jagged2DArray<T>::~Jagged2DArray()
	{
	}

	template <class T>
	void Jagged2DArray<T>::init_rectangular(const size_t nbr_rows, const size_t nbr_cols)
	{
		size_ = nbr_rows;
		row_starts_.reset(new size_t[nbr_rows + 1]);
		data_.reset(new T[nbr_rows * nbr_cols]);
		// remember that row_starts_ has 1 element more than we have rows
		for (size_t i = 0; i <= size_; ++i) {
			row_starts_[i] = i * nbr_cols;
		}
	}

	template <class T>
	void Jagged2DArray<T>::init_triangular(const size_t nbr_rows, const bool include_diagonal)
	{
		size_ = nbr_rows;
		row_starts_.reset(new size_t[nbr_rows + 1]);
		data_.reset(new T[triangular_size(nbr_rows, include_diagonal)]);
		size_t start = 0;
		for (size_t i = 0; i <= size_; ++i) {
			row_starts_[i] = start;
			start += include_diagonal ? (i + 1) : i;
		}
	}

	template <class T>
	size_t Jagged2DArray<T>::triangular_size(size_t N, bool include_diagonal)
	{
		return include_diagonal ? (N*(N+1))/2 : (N*(N-1))/2;
	}

	template <class T>
	Jagged2DArray<T>& Jagged2DArray<T>::operator=(const Jagged2DArray<T>& other)
	{
		// check for self-assingment
		if (this != &other) {
			// use copy-and-swap technique
			Jagged2DArray<T> copy(other);
			size_ = copy.size_;
			row_starts_.swap(copy.row_starts_);
			data_.swap(copy.data_);
		}
		return *this;
	}

	template <class T>
	void Jagged2DArray<T>::swap(Jagged2DArray<T>& other)
	{
		if (this != &other) {
			using std::swap;
			data_.swap(other.data_);
			row_starts_.swap(other.row_starts_);
			swap(size_, other.size_);
		}
	}

	template <class T>
	inline size_t Jagged2DArray<T>::size() const
	{
		return size_;
	}

	// Implementation detail
	inline std::string ja_err_msg_row_(const size_t row)
	{
		return std::string("Row index out of range: " + stringify(row));
	}

	template <class T>
	inline size_t Jagged2DArray<T>::row_size(const size_t row) const
	{
	#ifndef NDEBUG
		if (row >= size_) {
			throw std::out_of_range(ja_err_msg_row_(row));
		}
	#endif
		return static_cast<size_t>(row_starts_[row + 1] - row_starts_[row]);
	}

	template <class T>
	inline size_t Jagged2DArray<T>::nbr_elements() const
	{
		return row_starts_[size_];
	}

	template <class T>
	inline T& Jagged2DArray<T>::operator()(const size_t row, const size_t row_elem)
	{
	#ifndef NDEBUG
		if (row >= size_) {
			throw std::out_of_range(ja_err_msg_row_(row));
		}
		if (row_elem >= row_size(row)) {
			throw std::out_of_range("Row element index out of range");
		}
	#endif
		return data_[row_starts_[row] + row_elem];
	}

	template <class T>
	inline const T& Jagged2DArray<T>::operator()(const size_t row, const size_t row_elem) const
	{
	#ifndef NDEBUG
		if (row >= size_) {
			throw std::out_of_range(ja_err_msg_row_(row));
		}
		if (row_elem >= row_size(row)) {
			throw std::out_of_range("Row element index out of range");
		}
	#endif
		return data_[row_starts_[row] + row_elem];
	}

	template <class T>
	inline Jagged2DArrayRowRef<T> Jagged2DArray<T>::operator[](const size_t row_idx)
	{
	#ifndef NDEBUG
		if (row_idx >= size_) {
			throw std::out_of_range(ja_err_msg_row_(row_idx));
		}
	#endif
		return Jagged2DArrayRowRef<T>(data_.get() + row_starts_[row_idx], row_size(row_idx));
	}

	template <class T>
	inline const Jagged2DArrayRowRef<T> Jagged2DArray<T>::operator[](const size_t row_idx) const
	{
	#ifndef NDEBUG
		if (row_idx >= size_) {
			throw std::out_of_range(ja_err_msg_row_(row_idx));
		}
	#endif
		return Jagged2DArrayRowRef<T>(data_.get() + row_starts_[row_idx], row_size(row_idx));
	}

	template <class T>
	void Jagged2DArray<T>::fill(const T& elem)
	{
		std::fill(data_.get(), data_.get() + nbr_elements(), elem);
	}

	template <class T>
	inline Jagged2DArrayRowRef<T> Jagged2DArray<T>::flat_form()
	{
		return Jagged2DArrayRowRef<T>(data_.get(), nbr_elements());
	}

	template <class T>
	inline const Jagged2DArrayRowRef<T> Jagged2DArray<T>::flat_form() const
	{
		return Jagged2DArrayRowRef<T>(data_.get(), nbr_elements());
	}

	template <class T>
	inline typename Jagged2DArrayRowRef<T>::iterator Jagged2DArray<T>::flat_begin()
	{
		return data_.get();
	}

	template <class T>
	inline typename Jagged2DArrayRowRef<T>::iterator Jagged2DArray<T>::flat_end()
	{
		return data_.get() + nbr_elements();
	}

	template <class T>
	inline typename Jagged2DArrayRowRef<T>::const_iterator Jagged2DArray<T>::flat_begin() const
	{
		return data_.get();
	}

	template <class T>
	inline typename Jagged2DArrayRowRef<T>::const_iterator Jagged2DArray<T>::flat_end() const
	{
		return data_.get() + nbr_elements();
	}

	// Free functions for jagged arrays

	/** Does not throw */
	template <class T>
	void swap(Jagged2DArray<T>& a, Jagged2DArray<T>& b)
	{
		a.swap(b);
	}

	/** If it fails, it does not change the input object */
	template <class T>
	void resize(Jagged2DArray<T>& a, size_t nbr_rows, size_t nbr_cols)
	{
		Jagged2DArray<T> nowy(nbr_rows, nbr_cols);
		a.swap(nowy);
	}

	/** If it fails, it does not change the input object. Caller must provide iterators begin() and end()
	with new row sizes.
	@tparam Iter const iterator
	*/
	template <class T, class Iter>
	void resize(Jagged2DArray<T>& a, Iter begin, Iter end)
	{
		Jagged2DArray<T> nowy(begin, end);
		a.swap(nowy);
	}

	template <class T>
	void print(const Jagged2DArray<T>& array, std::ostream& out)
	{
		out << "Nbr rows: " << array.size() << std::endl;
		for (size_t i = 0; i < array.size(); ++i) {
			assert( array.row_size(i) == array[i].size() );
			out << "[" << array.row_size(i) << "] ";
			for (size_t j = 0; j < array.row_size(i); ++j) {
				out << array(i, j) << " ";
			}
			out << std::endl;
		}
	}
}}
