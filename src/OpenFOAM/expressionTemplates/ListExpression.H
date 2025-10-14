/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024,2025 M. Janssens
    Copyright (C) 2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::ListExpression

Description
    Expression templates for List

SourceFiles
    ListExpression.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListExpression_H
#define Foam_ListExpression_H

#include <cassert>
#include <execution>
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Expression
{

/*---------------------------------------------------------------------------*\
                       Namespace Expression Declarations
\*---------------------------------------------------------------------------*/

template<typename E, typename ValType>
class ListExpression
{
public:
    static constexpr bool is_leaf = false;

    ValType operator[](const label i) const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this)[i];
    }
    auto size() const noexcept { return static_cast<E const&>(*this).size(); }

    //- Helper: assign to passed in list
    template<class Container>
    Container& evaluate(Container& lst) const
    {
        lst.resize_nocopy(this->size());

#ifdef _OPENMP
        int loopLen = this->size();
        _Pragma("omp target teams distribute parallel for if(loopLen > 1000)")
        for(int i = 0; i < loopLen; ++i)
        {
            auto src = static_cast<E const&>(*this).cbegin() + i;
            lst[i] = *src;
        }
#else
    std::copy
    (
        #if defined(__APPLE__) && defined(__clang__)
        // Currently not working (2025-10-06) with par_unseq. Works with
        // sequential.
        #else
        std::execution::par_unseq,
        #endif
        static_cast<E const&>(*this).cbegin(),
        static_cast<E const&>(*this).cend(),
        lst.begin()
    );
#endif
        return lst;
    }
};


// Wrap of List
// ~~~~~~~~~~~~

template<class T>
class ListWrap
:
    public ListExpression<ListWrap<T>, T>
{
    List<T> elems_;

public:
    static constexpr bool is_leaf = true;

    typedef T value_type;

    //- Copy construct
    ListWrap(const ListWrap<T>& elems)
    :
        elems_(elems.data())
    {}

    //- Move construct
    ListWrap(ListWrap<T>&& elems)
    :
        elems_(std::move(elems.data()))
    {}

    // Construct from initializer list 
    ListWrap(std::initializer_list<T> elems)
    :
        elems_(elems)
    {}

    ListWrap(const UList<T>& elems)
    :
        elems_(elems)
    {}

    // Construct from ListExpression, forcing its evaluation.
    template<typename E>
    ListWrap(const ListExpression<E, typename E::value_type>& expr)
    :
        elems_(expr.size())
    {
        expr.evaluate(elems_);
    }

    // Access underlying data
    const auto& data() const
    { 
        return elems_; 
    }

    // Access underlying data
    auto& data()
    { 
        return elems_; 
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    value_type& operator[](const label i)
    {
        return elems_[i];
    }

    auto size() const noexcept
    {
        return elems_.size();
    }

    typedef typename List<T>::const_iterator const_iterator;
    auto cbegin() const
    {
        return elems_.cbegin();
    }
    auto cend() const
    {
        return elems_.cend();
    }

    typedef typename List<T>::iterator iterator;
    auto begin()
    {
        return elems_.begin();
    }
    auto end()
    {
        return elems_.end();
    }
};


// Wrap of non-const reference to List
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<class T>
class ListRefWrap
:
    public ListExpression<ListRefWrap<T>, T>
{
    List<T>& elems_;

public:
    static constexpr bool is_leaf = false;  //true;

    typedef T value_type;

    //- Copy construct
    ListRefWrap(ListRefWrap<T>& w)
    :
        elems_(w.elems_)
    {}

    //- Move construct
    ListRefWrap(ListRefWrap<T>&& elems)
    :
        elems_(elems)
    {}

    //- Construct from List and expected size. This is to enforce that
    //- our size() member function returns the correct size. Use the
    //- construct-from-expression otherwise.
    ListRefWrap(const label size, List<T>& elems)
    :
        elems_(elems)
    {
        elems_.resize_nocopy(size);
    }

    // Construct from ListExpression, forcing its evaluation.
    template<typename E>
    ListRefWrap
    (
        List<T>& elems,
        const ListExpression<E, typename E::value_type>& expr
    )
    :
        elems_(elems)
    {
        expr.evaluate(elems_);
    }

    //- Assignment
    template<typename E>
    void operator=(const ListExpression<E, typename E::value_type>& expr)
    {
        expr.evaluate(elems_);
    }

    // Access underlying data
    const auto& data() const
    { 
        return elems_; 
    }

    // Access underlying data
    auto& data()
    { 
        return elems_; 
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    value_type& operator[](const label i)
    {
        return elems_[i];
    }

    auto size() const noexcept
    {
        return elems_.size();
    }

    typedef typename List<T>::const_iterator const_iterator;
    auto cbegin() const
    {
        return elems_.cbegin();
    }
    auto cend() const
    {
        return elems_.cend();
    }

    typedef typename List<T>::iterator iterator;
    auto begin()
    {
        return elems_.begin();
    }
    auto end()
    {
        return elems_.end();
    }
};


// Wrap of const reference to List
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<class T>
class ListConstRefWrap
:
    public ListExpression<ListConstRefWrap<T>, T>
{
    const List<T>& elems_;

public:

    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    typedef T value_type;

    // construct from components
    ListConstRefWrap(const List<T>& elems)
    :
        elems_(elems)
    {}

    // returns the underlying data
    const auto& data() const
    { 
        return elems_; 
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    auto size() const noexcept
    {
        return elems_.size();
    }

    //typedef typename List<T>::const_iterator const_iterator;
    struct const_iterator
    {
        // Minimalistic const_iter implementation
        typedef typename List<T>::const_iterator ConstIter;
        typedef typename std::iterator_traits<ConstIter> ConstIterTraits;
        using difference_type = typename ConstIterTraits::difference_type;
        using value_type = typename ConstIterTraits::value_type;
        using pointer = typename ConstIterTraits::pointer;
        using reference = typename ConstIterTraits::reference;
//using iterator_category = typename ConstIterTraits::iterator_category;
        using iterator_category = std::random_access_iterator_tag;

        ConstIter base_;

        const_iterator(ConstIter it) : base_(it) {}
        const_iterator() : base_() {}
        ConstIter& operator=(ConstIter it)
        {
            base_ = it;
            return base_;
        }

    // Indexing operators

        value_type operator*() const
        {
            return *base_;
        }
        //value_type operator->() const
        //{
        //    return *base_;
        //}
        value_type operator[](const difference_type i) const
        {
           return *(base_+i);
        }

    // Increment/Decrement

        const_iterator& operator++()
        {
            ++base_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++base_;
            return old;
        }
        const_iterator& operator--()
        {
            --base_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --base_;
            return old;
        }

    // Arithmetic operators

        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(base_+i);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(base_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (base_ - it.base_);
        }

    // Comparison operators

        bool operator==(const const_iterator& it) const
        {
            return (base_ == it.base_);
        }
        bool operator!=(const const_iterator& it) const
        {
            return (base_ != it.base_);
        }
        bool operator<(const const_iterator& it) const
        {
            return (base_ < it.base_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (base_ <= it.base_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (base_ > it.base_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (base_ >= it.base_);
        }

    // Compound assignment 

        const_iterator& operator+=(const difference_type n)
        {
            base_ += n;
            return *this;
        }
        const_iterator& operator-=(const difference_type n)
        {
            base_ -= n;
            return *this;
        }
    };
    auto cbegin() const
    {
        return const_iterator(elems_.cbegin());
    }
    auto cend() const
    {
        return const_iterator(elems_.cend());
    }
};


// Wrap of tmp to List
// ~~~~~~~~~~~~~~~~~~~

template<class T>
class ListTmpWrap
:
    public ListExpression<ListTmpWrap<T>, typename T::value_type>
{
    tmp<T> elems_;

public:

    //- Have expressions use copy to maintain tmp refCount
    static constexpr bool is_leaf = false;

    typedef typename T::value_type value_type;

    // Construct from components
    ListTmpWrap(const tmp<T>& elems)
    :
        elems_(elems)
    {}

    // Construct from components
    ListTmpWrap(tmp<T>& elems)
    :
        elems_(elems)
    {}

    // Construct from ListExpression, forcing its evaluation.
    template<typename E>
    ListTmpWrap
    (
        tmp<T>& elems,
        const ListExpression<E, typename E::value_type>& expr
    )
    :
        elems_(elems)
    {
        expr.evaluate(elems_.ref());
    }

    //- Assignment
    template<typename E>
    void operator=(const ListExpression<E, typename E::value_type>& expr)
    {
        expr.evaluate(elems_.ref());
    }

    // returns the underlying data
    const T& data() const
    {
        return elems_();
    }

    value_type operator[](const label i) const
    {
        return elems_()[i];
    }

    auto size() const noexcept
    {
        return elems_().size();
    }

    typedef typename T::const_iterator const_iterator;
    auto cbegin() const
    {
        return elems_().cbegin();
    }
    auto cend() const
    {
        return elems_().cend();
    }
    typedef typename List<T>::iterator iterator;
    auto begin()
    {
        return elems_().begin();
    }
    auto end()
    {
        return elems_().end();
    }
};


template<class T>
class ListConstTmpWrap
:
    public ListExpression<ListConstTmpWrap<T>, typename T::value_type>
{
    const tmp<T> elems_;

public:

    //- Have expressions use copy to maintain tmp refCount
    static constexpr bool is_leaf = false;

    typedef typename T::value_type value_type;

    // Construct from components
    ListConstTmpWrap(const tmp<T>& elems)
    :
        elems_(elems)
    {}

    // returns the underlying data
    const T& data() const
    {
        return elems_();
    }

    value_type operator[](const label i) const
    {
        return elems_()[i];
    }

    auto size() const noexcept
    {
        return elems_().size();
    }

    //- Define our own iterator so we can refer to it. 'normal', e.g. double*
    //  iterator will not work we'd have to use the iterator traits to get
    //  at its e.g. difference_type from other classes. TBD.
    struct const_iterator
    {
        // Minimalistic const_iter implementation
        typedef typename T::const_iterator ConstIter;
        typedef typename std::iterator_traits<ConstIter> ConstIterTraits;
        using difference_type = typename ConstIterTraits::difference_type;
        using value_type = typename ConstIterTraits::value_type;
        using pointer = typename ConstIterTraits::pointer;
        using reference = typename ConstIterTraits::reference;
        using iterator_category = typename ConstIterTraits::iterator_category;

        ConstIter base_;

        const_iterator(ConstIter it) : base_(it) {}
        value_type operator*() const
        {
            return *base_;
        }
        value_type operator[](const difference_type i) const
        {
            return base_[i];
        }
        const_iterator& operator++()
        {
            ++base_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++base_;
            return old;
        }
        const_iterator& operator--()
        {
            --base_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --base_;
            return old;
        }
        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(base_+i);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(base_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (base_ - it.base_);
        }
        bool operator==(const const_iterator& it) const
        {
            return (base_ == it.base_);
        }
        bool operator!=(const const_iterator& it) const
        {
            return (base_ != it.base_);
        }
        bool operator<(const const_iterator& it) const
        {
            return (base_ < it.base_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (base_ <= it.base_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (base_ > it.base_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (base_ >= it.base_);
        }
        const_iterator& operator+=(const difference_type n)
        {
            base_ += n;
            return base_;
        }
        const_iterator& operator-=(const difference_type n)
        {
            base_ -= n;
            return base_;
        }
    };
    //typedef typename T::const_iterator const_iterator;
    auto cbegin() const
    {
        return const_iterator(elems_().cbegin());
    }
    auto cend() const
    {
        return const_iterator(elems_().cend());
    }
};


// Wrap of indirection
// ~~~~~~~~~~~~~~~~~~~

// E1 : expression to get value of element i
template<typename E1>
class IndirectConstWrap
:
    public ListExpression
    <
        IndirectConstWrap<E1>,
        typename E1::value_type
    >
{
    //- Expression for cell values
    typename std::conditional
    <
        E1::is_leaf,
        const E1&,
        const E1
    >::type values_;

    const labelUList& addr_;

public:

    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    typedef typename E1::value_type value_type;

    // construct from components
    IndirectConstWrap
    (
        const E1& values,
        const labelUList& addr
    )
    :
        values_(values),
        addr_(addr)
    {}

    value_type operator[](const label i) const
    {
        return values_[addr_[i]];
    }
    label size() const
    {
        return addr_.size();
    }

    struct const_iterator
    {
        // See IndirectListBase::const_iterator

        //- Step type for addressing
        using difference_type = label;

        //- Iterator definition for E1
        typedef typename E1::const_iterator ConstIter;
        //using difference_type = typename ConstIter::difference_type;
        using value_type = typename ConstIter::value_type;
        using pointer = typename ConstIter::pointer;
        using reference = typename ConstIter::reference;
        using iterator_category = std::random_access_iterator_tag;

        ConstIter begin_;
        typename labelUList::const_iterator iter_;

        const_iterator
        (
            ConstIter begin,
            typename labelUList::const_iterator addrIter
        )
        :
            begin_(begin),
            iter_(addrIter)
        {}
        auto operator*() const
        {
            return *(begin_ + *iter_);
        }
        value_type operator[](const difference_type i) const
        {
            return *(begin_+iter_[i]);
        }
        const_iterator& operator++()
        {
            ++iter_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++iter_;
            return old;
        }
        const_iterator& operator--()
        {
            --iter_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --iter_;
            return old;
        }
        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(begin_, iter_+i);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(begin_, iter_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (iter_ - it.iter_);
        }
        bool operator==(const const_iterator& it) const
        {
            return iter_ == it.iter_;
        }
        bool operator!=(const const_iterator& it) const
        {
            return iter_ != it.iter_;
        }
        bool operator<(const const_iterator& it) const
        {
            return (iter_ < it.iter_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (iter_ <= it.iter_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (iter_ > it.iter_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (iter_ >= it.iter_);
        }
        const_iterator& operator+=(const difference_type n)
        {
            iter_ += n;
            return *this;
        }
        const_iterator& operator-=(const difference_type n)
        {
            iter_ -= n;
            return *this;
        }
    };
    auto cbegin() const
    {
        return const_iterator(values_.cbegin(), addr_.cbegin());
    }
    auto cend() const
    {
        return const_iterator(values_.cbegin(), addr_.cend());
    }
};


// Expresssions on List
// ~~~~~~~~~~~~~~~~~~~~

#define EXPRESSION_FUNCTION1(Func, BaseFunc, OpFunc)                           \
template<typename E1>                                                          \
class OpFunc                                                                   \
:                                                                              \
    public ListExpression<OpFunc<E1>, typename E1::value_type>                 \
{                                                                              \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    OpFunc(const E1& u)                                                        \
    :                                                                          \
        u_(u)                                                                  \
    {}                                                                         \
                                                                               \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return BaseFunc(u_[i]);                                                \
    }                                                                          \
                                                                               \
    auto size() const noexcept { return u_.size(); }                           \
                                                                               \
    struct const_iterator                                                      \
    {                                                                          \
        typedef typename E1::const_iterator ConstIter;                         \
        using difference_type = typename ConstIter::difference_type;           \
        using value_type = typename ConstIter::value_type;                     \
        using pointer = typename ConstIter::pointer;                           \
        using reference = typename ConstIter::reference;                       \
        using iterator_category = typename ConstIter::iterator_category;       \
                                                                               \
        ConstIter base_;                                                       \
                                                                               \
        const_iterator(ConstIter it) : base_(it) {}                            \
        auto operator*() const                                           \
        {                                                                      \
            return BaseFunc(*base_);                                           \
        }                                                                      \
        auto operator[](const difference_type i) const                   \
        {                                                                      \
            return BaseFunc(base_[i]);                                         \
        }                                                                      \
        const_iterator& operator++()                                           \
        {                                                                      \
            ++base_;                                                           \
            return *this;                                                      \
        }                                                                      \
        const_iterator operator++(int)                                         \
        {                                                                      \
            const_iterator old(*this);                                         \
            ++base_;                                                           \
            return old;                                                        \
        }                                                                      \
        const_iterator& operator--()                                           \
        {                                                                      \
            --base_;                                                           \
            return *this;                                                      \
        }                                                                      \
        const_iterator operator--(int)                                         \
        {                                                                      \
            const_iterator old(*this);                                         \
            --base_;                                                           \
            return old;                                                        \
        }                                                                      \
        const_iterator operator+(const difference_type i) const                \
        {                                                                      \
            return const_iterator(base_+i);                                    \
        }                                                                      \
        const_iterator operator-(const difference_type i) const                \
        {                                                                      \
            return const_iterator(base_-i);                                    \
        }                                                                      \
        difference_type operator-(const const_iterator& it) const              \
        {                                                                      \
            return (base_ - it.base_);                                         \
        }                                                                      \
        bool operator==(const const_iterator& it) const                        \
        {                                                                      \
            return base_ == it.base_;                                          \
        }                                                                      \
        bool operator!=(const const_iterator& it) const                        \
        {                                                                      \
            return base_ != it.base_;                                          \
        }                                                                      \
        bool operator<(const const_iterator& it) const                         \
        {                                                                      \
            return (base_ < it.base_);                                         \
        }                                                                      \
        bool operator<=(const const_iterator& it) const                        \
        {                                                                      \
            return (base_ <= it.base_);                                        \
        }                                                                      \
        bool operator>(const const_iterator& it) const                         \
        {                                                                      \
            return (base_ > it.base_);                                         \
        }                                                                      \
        bool operator>=(const const_iterator& it) const                        \
        {                                                                      \
            return (base_ >= it.base_);                                        \
        }                                                                      \
        const_iterator& operator+=(const difference_type n)                    \
        {                                                                      \
            base_ += n;                                                        \
            return *this;                                                      \
        }                                                                      \
        const_iterator& operator-=(const difference_type n)                    \
        {                                                                      \
            base_ -= n;                                                        \
            return *this;                                                      \
        }                                                                      \
    };                                                                         \
    const_iterator cbegin() const                                              \
    {                                                                          \
        return const_iterator(u_.cbegin());                                    \
    }                                                                          \
    const_iterator cend() const                                                \
    {                                                                          \
        return const_iterator(u_.cend());                                      \
    }                                                                          \
};                                                                             \
template<typename E1>                                                          \
OpFunc<E1> Func                                                                \
(                                                                              \
    const ListExpression<E1, typename E1::value_type>& u                       \
)                                                                              \
{                                                                              \
    return OpFunc<E1>(static_cast<const E1&>(u));                              \
}


// Do '-' separately until we work out macro expansion...
template<typename E1>
class List_negate
:
    public ListExpression<List_negate<E1>, typename E1::value_type>
{
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

public:
    static constexpr bool is_leaf = false;

    typedef typename E1::value_type value_type;

    List_negate(const E1& u)
    :
        u_(u)
    {}

    auto operator[](const label i) const
    {
        return -u_[i];
    }

    auto size() const noexcept { return u_.size(); }

    struct const_iterator
    {
        // Minimalistic const_iter implementation

        typedef typename E1::const_iterator ConstIter;
        using difference_type = typename ConstIter::difference_type;
        using value_type = typename ConstIter::value_type;
        using pointer = typename ConstIter::pointer;
        using reference = typename ConstIter::reference;
        using iterator_category = typename ConstIter::iterator_category;

        ConstIter base_;

        const_iterator(ConstIter it) : base_(it) {}
        value_type operator*() const
        {
            return -(*base_);
        }
        value_type operator[](const difference_type i) const
        {
           return -*(base_+i);
        }
        const_iterator& operator++()
        {
            ++base_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++base_;
            return old;
        }
        const_iterator& operator--()
        {
            --base_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --base_;
            return old;
        }
        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(base_+i);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(base_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (base_ - it.base_);
        }
        bool operator==(const const_iterator& it) const
        {
            return base_ == it.base_;
        }
        bool operator!=(const const_iterator& it) const
        {
            return base_ != it.base_;
        }
        bool operator<(const const_iterator& it) const
        {
            return (base_ < it.base_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (base_ <= it.base_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (base_ > it.base_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (base_ >= it.base_);
        }
        const_iterator& operator+=(const difference_type n)
        {
            base_ += n;
            return base_;
        }
        const_iterator& operator-=(const difference_type n)
        {
            base_ -= n;
            return base_;
        }
    };
    const_iterator cbegin() const
    {
        return const_iterator(u_.cbegin());
    }
    const_iterator cend() const
    {
        return const_iterator(u_.cend());
    }
};
template<typename E1>
List_negate<E1> operator-
(
    const ListExpression<E1, typename E1::value_type>& u
)
{
    return List_negate<E1>(static_cast<const E1&>(u));
}


#define EXPRESSION_FUNCTION2(Func, BaseFunc, OpFunc)                           \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public ListExpression                                                      \
    <                                                                          \
        OpFunc<E1, E2>,                                                        \
        typename E1::value_type                                                \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    /* Type to return for internal field */                                    \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
     u_(u), v_(v)                                                              \
    {                                                                          \
        assert(u.size() == v.size());                                          \
    }                                                                          \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return BaseFunc(u_[i], v_[i]);                                         \
    }                                                                          \
    auto size() const noexcept {return u_.size(); }                            \
                                                                               \
    struct const_iterator                                                      \
    {                                                                          \
        typedef typename E1::const_iterator E1Iter;                            \
        typedef typename E2::const_iterator E2Iter;                            \
        using difference_type = typename E1Iter::difference_type;              \
        using value_type = typename E1Iter::value_type;                        \
        using pointer = typename E1Iter::pointer;                              \
        using reference = typename E1Iter::reference;                          \
        using iterator_category = typename E1Iter::iterator_category;          \
                                                                               \
        E1Iter uIter_;                                                         \
        E2Iter vIter_;                                                         \
                                                                               \
        const_iterator(E1Iter uIter, E2Iter vIter)                             \
        :                                                                      \
            uIter_(uIter),                                                     \
            vIter_(vIter)                                                      \
        {}                                                                     \
        auto operator*() const                                                   \
        {                                                                      \
            return BaseFunc(*uIter_, *vIter_);                                 \
        }                                                                      \
        const_iterator& operator++()                                           \
        {                                                                      \
            ++uIter_;                                                          \
            ++vIter_;                                                          \
            return *this;                                                      \
        }                                                                      \
        difference_type operator-(const const_iterator& it) const              \
        {                                                                      \
            return (uIter_ - it.uIter_);                                       \
        }                                                                      \
        const_iterator operator+(const label i) const                          \
        {                                                                      \
            return const_iterator(uIter_+i, vIter_+i);                         \
        }                                                                      \
        bool operator==(const const_iterator& it) const                        \
        {                                                                      \
            return uIter_ == it.uIter_ && vIter_ == it.vIter_;                 \
        }                                                                      \
        bool operator!=(const const_iterator& it) const                        \
        {                                                                      \
            return !operator==(it);                                            \
        }                                                                      \
        bool operator<(const const_iterator& it) const                         \
        {                                                                      \
            return (uIter_ < it.uIter_);                                       \
        }                                                                      \
    };                                                                         \
    const_iterator cbegin() const                                              \
    {                                                                          \
        return const_iterator(u_.cbegin(), v_.cbegin());                       \
    }                                                                          \
    const_iterator cend() const                                                \
    {                                                                          \
        return const_iterator(u_.cend(), v_.cend());                           \
    }                                                                          \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2> Func                                                            \
(                                                                              \
    ListExpression                                                             \
    <                                                                          \
        E1,                                                                    \
        typename E1::value_type                                                \
    > const& u,                                                                \
    ListExpression                                                             \
    <                                                                          \
        E2,                                                                    \
        typename E2::value_type                                                \
    > const& v                                                                 \
)                                                                              \
{                                                                              \
    return OpFunc<E1, E2>                                                      \
    (                                                                          \
        static_cast<const E1&>(u),                                             \
        static_cast<const E2&>(v)                                              \
    );                                                                         \
}
#define EXPRESSION_OPERATOR(Op, OpFunc)                                        \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public ListExpression                                                      \
    <                                                                          \
        OpFunc<E1, E2>,                                                        \
        typename E1::value_type                                                \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    /* Type to return for internal field */                                    \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
     u_(u), v_(v)                                                              \
    {                                                                          \
        assert(u.size() == v.size());                                          \
    }                                                                          \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return u_[i] Op v_[i];                                                 \
    }                                                                          \
    auto size() const noexcept { return u_.size(); }                           \
                                                                               \
    struct const_iterator                                                      \
    {                                                                          \
        typedef typename E1::const_iterator E1Iter;                            \
        typedef typename E2::const_iterator E2Iter;                            \
        using difference_type = typename E1Iter::difference_type;              \
        using value_type = typename E1Iter::value_type;                        \
        using pointer = typename E1Iter::pointer;                              \
        using reference = typename E1Iter::reference;                          \
        using iterator_category = typename E1Iter::iterator_category;          \
                                                                               \
        E1Iter uIter_;                                                         \
        E2Iter vIter_;                                                         \
                                                                               \
        const_iterator(E1Iter uIter, E2Iter vIter)                             \
        :                                                                      \
            uIter_(uIter),                                                     \
            vIter_(vIter)                                                      \
        {}                                                                     \
        value_type operator*() const                                           \
        {                                                                      \
            return *uIter_ Op *vIter_;                                         \
        }                                                                      \
        value_type operator[](const difference_type i) const                   \
        {                                                                      \
           return uIter_[i] Op vIter_[i];                                      \
        }                                                                      \
        const_iterator& operator++()                                           \
        {                                                                      \
            ++uIter_;                                                          \
            ++vIter_;                                                          \
            return *this;                                                      \
        }                                                                      \
        const_iterator operator++(int)                                         \
        {                                                                      \
            const_iterator old(*this);                                         \
            ++uIter_;                                                          \
            ++vIter_;                                                          \
            return old;                                                        \
        }                                                                      \
        const_iterator& operator--()                                           \
        {                                                                      \
            --uIter_;                                                          \
            --vIter_;                                                          \
            return *this;                                                      \
        }                                                                      \
        const_iterator operator--(int)                                         \
        {                                                                      \
            const_iterator old(*this);                                         \
            --uIter_;                                                          \
            --vIter_;                                                          \
            return old;                                                        \
        }                                                                      \
        const_iterator operator+(const difference_type i) const                \
        {                                                                      \
            return const_iterator(uIter_+i, vIter_+i);                         \
        }                                                                      \
        const_iterator operator-(const difference_type i) const                \
        {                                                                      \
            return const_iterator(uIter_-i, vIter_-i);                         \
        }                                                                      \
        difference_type operator-(const const_iterator& it) const              \
        {                                                                      \
            return (uIter_ - it.uIter_);                                       \
        }                                                                      \
        bool operator==(const const_iterator& it) const                        \
        {                                                                      \
            return uIter_ == it.uIter_ && vIter_ == it.vIter_;                 \
        }                                                                      \
        bool operator!=(const const_iterator& it) const                        \
        {                                                                      \
            return !operator==(it);                                            \
        }                                                                      \
        bool operator<(const const_iterator& it) const                         \
        {                                                                      \
            return (uIter_ < it.uIter_);                                       \
        }                                                                      \
        bool operator<=(const const_iterator& it) const                        \
        {                                                                      \
            return (uIter_ <= it.uIter_);                                      \
        }                                                                      \
        bool operator>(const const_iterator& it) const                         \
        {                                                                      \
            return (uIter_ > it.uIter_);                                       \
        }                                                                      \
        bool operator>=(const const_iterator& it) const                        \
        {                                                                      \
            return (uIter_ >= it.uIter_);                                      \
        }                                                                      \
        const_iterator& operator+=(const difference_type n)                    \
        {                                                                      \
            uIter_ += n;                                                       \
            vIter_ += n;                                                       \
            return *this;                                                      \
        }                                                                      \
        const_iterator& operator-=(const difference_type n)                    \
        {                                                                      \
            uIter_ -= n;                                                       \
            vIter_ -= n;                                                       \
            return *this;                                                      \
        }                                                                      \
    };                                                                         \
    const_iterator cbegin() const                                              \
    {                                                                          \
        return const_iterator(u_.cbegin(), v_.cbegin());                       \
    }                                                                          \
    const_iterator cend() const                                                \
    {                                                                          \
        return const_iterator(u_.cend(), v_.cend());                           \
    }                                                                          \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2>                                                                 \
operator Op                                                                    \
(                                                                              \
    ListExpression                                                             \
    <                                                                          \
        E1,                                                                    \
        typename E1::value_type                                                \
    > const& u,                                                                \
    ListExpression                                                             \
    <                                                                          \
        E2,                                                                    \
        typename E2::value_type                                                \
    > const& v                                                                 \
)                                                                              \
{                                                                              \
    return OpFunc<E1, E2>                                                      \
    (                                                                          \
        static_cast<const E1&>(u),                                             \
        static_cast<const E2&>(v)                                              \
    );                                                                         \
}


//- List functions
// macro arguments:
//  - 'sin' : name of new function (in ListExpression namespace)
//  - '::sin' : per-element function to call
//  - 'List_sin' : generated expression helper class
EXPRESSION_FUNCTION1(sin, ::sin, List_sin)
EXPRESSION_FUNCTION1(cos, ::cos, List_cos)
EXPRESSION_FUNCTION1(tan, ::tan, List_tan)
EXPRESSION_FUNCTION1(sinh, ::sinh, List_sinh)
EXPRESSION_FUNCTION1(cosh, ::cosh, List_cosh)
EXPRESSION_FUNCTION1(tanh, ::tanh, List_tanh)
EXPRESSION_FUNCTION1(sqr, Foam::sqr, List_sqr)
EXPRESSION_FUNCTION1(sqrt, Foam::sqrt, List_sqrt)
EXPRESSION_FUNCTION1(magSqr, Foam::magSqr, List_magSqr)
EXPRESSION_FUNCTION1(mag, Foam::mag, List_mag)
EXPRESSION_FUNCTION1(symm, Foam::symm, List_symm)
//EXPRESSION_FUNCTION1(-, -, List_negate)
EXPRESSION_FUNCTION1(pow2, Foam::pow2, List_pow2)
EXPRESSION_FUNCTION1(pow3, Foam::pow3, List_pow3)
EXPRESSION_FUNCTION1(pow4, Foam::pow4, List_pow4)

#undef EXPRESSION_FUNCTION1

EXPRESSION_FUNCTION2(min, Foam::min, List_min)
EXPRESSION_FUNCTION2(max, Foam::max, List_max)

#undef EXPRESSION_FUNCTION2

//- List operators
EXPRESSION_OPERATOR(+, List_add)
EXPRESSION_OPERATOR(-, List_subtract)
EXPRESSION_OPERATOR(*, List_multiply)
EXPRESSION_OPERATOR(/, List_divide)
EXPRESSION_OPERATOR(&, List_dot)

#undef EXPRESSION_OPERATOR

//- Not possible: a *= b has to return a new type so it can never be in-place
//EXPRESSION_LIST_OPERATOR(+=, List_inplace_add)
//EXPRESSION_LIST_OPERATOR(-=, List_inplace_subtract)
//EXPRESSION_LIST_OPERATOR(*=, List_inplace_multiply)
//EXPRESSION_LIST_OPERATOR(/=, List_inplace_divide)


// Interpolation
// ~~~~~~~~~~~~~

// E1 : expression to get value of element i
// E2 : expression to get value of element i
// E3 : expression to get interpolation weight
template<typename E1, typename E2, typename E3>
class List_interpolate
:
    public ListExpression
    <
        List_interpolate<E1, E2, E3>,
        typename E1::value_type
    >
{
    //- Expression for cell values1
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

    //- Expression for cell values2
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

    //- Expression for interpolation weight
    typename std::conditional<E3::is_leaf, const E3&, const E3>::type weight_;


public:

    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    typedef typename E1::value_type value_type;

    // construct from components
    List_interpolate(const E1& u, const E2& v, const E3& weight)
    :
        u_(u),
        v_(v),
        weight_(weight)
    {}

    value_type operator[](const label i) const
    {
        const auto ownVal = u_[i];
        const auto neiVal = v_[i];
        return weight_[i]*(ownVal-neiVal) + neiVal;
    }
    label size() const
    {
        return weight_.size();
    }

    struct const_iterator
    {
        typedef typename E1::const_iterator E1Iter;
        typedef typename E2::const_iterator E2Iter;
        typedef typename E3::const_iterator E3Iter;
        using difference_type = typename E3Iter::difference_type;
        using value_type = typename E1Iter::value_type;
        using pointer = typename E1Iter::pointer;
        using reference = typename E1Iter::reference;
        using iterator_category = typename E3Iter::iterator_category;
        E1Iter uIter_;
        E2Iter vIter_;
        E3Iter weightIter_;

        const_iterator(E1Iter uIter, E2Iter vIter, E3Iter weightIter)
        :
            uIter_(uIter),
            vIter_(vIter),
            weightIter_(weightIter)
        {}
        value_type operator*() const
        {
            const auto ownVal = *uIter_;
            const auto neiVal = *vIter_;
            return *weightIter_*(ownVal-neiVal) + neiVal;
        }
        value_type operator[](const difference_type i) const
        {
            const auto ownVal = uIter_[i];
            const auto neiVal = vIter_[i];
            return weightIter_[i]*(ownVal-neiVal) + neiVal;
        }
        const_iterator& operator++()
        {
            ++uIter_;
            ++vIter_;
            ++weightIter_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++uIter_;
            ++vIter_;
            ++weightIter_;
            return old;
        }
        const_iterator& operator--()
        {
            --uIter_;
            --vIter_;
            --weightIter_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --uIter_;
            --vIter_;
            --weightIter_;
            return old;
        }
        const_iterator operator+(const label i) const
        {
            return const_iterator(uIter_+i, vIter_+i, weightIter_+i);
        }
        const_iterator operator-(const label i) const
        {
            return const_iterator(uIter_-i, vIter_-i, weightIter_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (uIter_ - it.uIter_);
        }
        bool operator==(const const_iterator& it) const
        {
            return
                uIter_ == it.uIter_
             && vIter_ == it.vIter_
             && weightIter_ == it.weightIter_;
        }
        bool operator!=(const const_iterator& it) const
        {
            return !operator==(it);
        }
        bool operator<(const const_iterator& it) const
        {
            return (uIter_ < it.uIter_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (uIter_ <= it.uIter_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (uIter_ > it.uIter_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (uIter_ >= it.uIter_);
        }
        const_iterator& operator+=(const difference_type n)
        {
            uIter_ += n;
            vIter_ += n;
            weightIter_ += n;
            return *this;
        }
        const_iterator& operator-=(const difference_type n)
        {
            uIter_ -= n;
            vIter_ -= n;
            weightIter_ -= n;
            return *this;
        }
    };
    const_iterator cbegin() const
    {
        return const_iterator(u_.cbegin(), v_.cbegin(), weight_.cbegin());
    }
    const_iterator cend() const
    {
        return const_iterator(u_.cend(), v_.cend(), weight_.cend());
    }
};


// Expressions on constants
// ~~~~~~~~~~~~~~~~~~~~~~~~

template<class T>
class UniformListWrap
:
    public ListExpression<UniformListWrap<T>, T>
{
private:

    const label size_;

    const T val_;


public:

    static constexpr bool is_leaf = false;

    typedef T value_type;

    // construct from components
    UniformListWrap(const label size, const T val)
    :
        size_(size),
        val_(val)
    {}

    // construct from components
    UniformListWrap(const label size, const dimensioned<T>& val)
    :
        size_(size),
        val_(val.value())
    {}

    T operator[](const label i) const
    {
        return val_;
    }
    auto size() const noexcept
    {
        return size_;
    }

    struct const_iterator
    {
        // Minimalistic const_iter implementation
        using difference_type = label;
        using value_type = const T;
        using pointer = const T*;
        using reference = const T&;
        using iterator_category = std::random_access_iterator_tag;

        difference_type count_;

        const T val_;

        const_iterator(const difference_type count, const T val)
        :
            count_(count),
            val_(val)
        {}
        value_type operator*() const
        {
            return val_;
        }
        value_type operator[](const difference_type i) const
        {
           return val_;
        }
        const_iterator& operator++()
        {
            count_++;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++count_;
            return old;
        }
        const_iterator& operator--()
        {
            count_--;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --count_;
            return old;
        }
        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(count_+i, val_);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(count_-i, val_);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (count_ - it.count_);
        }
        bool operator==(const const_iterator& it) const
        {
            return val_ == it.val_ && count_ == it.count_;
        }
        bool operator!=(const const_iterator& it) const
        {
            return !operator==(it);
        }
        bool operator<(const const_iterator& it) const
        {
            return (count_ < it.count_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (count_ <= it.count_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (count_ > it.count_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (count_ >= it.count_);
        }
        const_iterator& operator+=(const difference_type n)
        {
            count_ += n;
            return *this;
        }
        const_iterator& operator-=(const difference_type n)
        {
            count_ -= n;
            return *this;
        }
    };
    auto cbegin() const
    {
        return const_iterator(0, val_);
    }
    auto cend() const
    {
        return const_iterator(size_, val_);
    }
};


//- Wrapper of multiple lists (e.g. FieldField)
template <class V>
class ListsConstRefWrap
:
    public ListExpression<ListsConstRefWrap<V>, typename V::value_type>
{
    List<const V*> data;
    List<size_t> offsets;
    size_t total_size = 0;

    void init()
    {
        offsets.resize_nocopy(data.size());
        total_size = 0;
        for (label i = 0; i < data.size(); ++i)
        {
            offsets[i] = total_size;
            total_size += data[i]->size();
        }
    }

public:
    static constexpr bool is_leaf = false;  //true;

    using value_type = typename V::value_type;  //T;
    typedef ListsConstRefWrap<V> this_type;

    ListsConstRefWrap() = delete;

    ListsConstRefWrap(const UList<const V*>& d)
    :
        data(d)
    {
        init();
    }
    ListsConstRefWrap(List<const V*>&& d)
    :
        data(d)
    {
        init();
    }
    template<class Container>
    ListsConstRefWrap(const Container& ds)
    {
        data.resize_nocopy(ds.size());
        forAll(ds, i)
        {
            const V& vals = static_cast<const V&>(ds[i]);
            data[i] = &vals;
        }
        init();
    }
    template<class Container>
    ListsConstRefWrap(const Container& ds, const labelUList& elems)
    {
        data.resize_nocopy(elems.size());
        forAll(ds, i)
        {
            const V& vals = static_cast<const V&>(ds[elems[i]]);
            data[i] = &vals;
        }
        init();
    }

    void push_back(const V* vec)
    {
        data.push_back(vec);
        offsets.push_back(total_size);
        total_size += vec->size();
    }

    value_type operator[](label i) const
    {    
        const label vec_index =
            std::upper_bound
            (
                offsets.begin(),
                offsets.end(),
                i
            )
         - offsets.begin() - 1;
        const label elem_index = i - offsets[vec_index];
        return (*data[vec_index])[elem_index];
    }
    label size() const { return total_size; }


    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef label difference_type;

    struct const_iterator
    {
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = this_type::difference_type;
        using value_type  = this_type::value_type;
        using pointer     = this_type::const_pointer;
        using reference   = this_type::const_reference;

        const this_type& base_; 
        label index_;

        const_iterator(const this_type& base, const label index)
        :
            base_(base),
            index_(index)
        {}
        const_iterator(const const_iterator& it)
        :
            base_(it.base_),
            index_(it.index_)
        {}
        const_iterator() = delete;  //: base_(0) {}
        const_iterator& operator=(const_iterator& it) = default;

    // Indexing operators

        value_type operator*() const
        {
            return base_[index_];
        }
        //value_type operator->() const
        //{
        //    return *base_;
        //}
        value_type operator[](const difference_type i) const
        {
           return base_[index_];
        }

    // Increment/Decrement

        const_iterator& operator++()
        {
            ++index_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++index_;
            return old;
        }
        const_iterator& operator--()
        {
            --index_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --index_;
            return old;
        }

    // Arithmetic operators

        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(base_, index_+i);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(base_, index_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (index_ - it.index_);
        }

    // Comparison operators

        bool operator==(const const_iterator& it) const
        {
            return (index_ == it.index_);
        }
        bool operator!=(const const_iterator& it) const
        {
            return (index_ != it.index_);
        }
        bool operator<(const const_iterator& it) const
        {
            return (index_ < it.index_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (index_ <= it.index_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (index_ > it.index_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (index_ >= it.index_);
        }

    // Compound assignment 

        const_iterator& operator+=(const difference_type n)
        {
            index_ += n;
            return *this;
        }
        const_iterator& operator-=(const difference_type n)
        {
            index_ -= n;
            return *this;
        }
    };
    auto cbegin() const
    {
        return const_iterator(*this, 0);
    }
    auto cend() const
    {
        return const_iterator(*this, this->size());
    }
};


template <class V>
class ListsRefWrap
:
    public ListExpression<ListsRefWrap<V>, typename V::value_type>//,T>, T>
{
    List<V*> data;
    List<size_t> offsets;
    size_t total_size = 0;

    void init()
    {
        offsets.resize_nocopy(data.size());
        total_size = 0;
        for (label i = 0; i < data.size(); ++i)
        {
            offsets[i] = total_size;
            total_size += data[i]->size();
        }
    }
    labelPair whichVecElem(const label i) const
    {
        const label vec_index =
        std::upper_bound
        (
            offsets.begin(),
            offsets.end(),
            i
        )
        - offsets.begin() - 1;
        const label elem_index = i - offsets[vec_index];
        return labelPair(vec_index, elem_index);
    }

public:
    static constexpr bool is_leaf = false;  //true;
    
    using value_type = typename V::value_type;  //T;
    typedef ListsRefWrap<V> this_type;

    ListsRefWrap() = delete;
    
    ListsRefWrap(UList<V*>& d)
    :
        data(d)
    {
        init();
    }
    ListsRefWrap(List<V*>&& d)
    :
        data(d)
    {
        init();
    }
    template<class Container>
    ListsRefWrap(Container& ds)
    {
        data.resize_nocopy(ds.size());
        forAll(ds, i)
        {
            V& vals = static_cast<V&>(ds[i]);
            data[i] = &vals;
        }
        init();
    }
    template<class Container>
    ListsRefWrap(Container& ds, const labelUList& elems)
    {
        data.resize_nocopy(elems.size());
        forAll(ds, i)
        {
            V& vals = static_cast<V&>(ds[i]);
            data[i] = &vals;
        }
        init();
    }

    void push_back(V* vec)
    {
        data.push_back(vec);
        offsets.push_back(total_size);
        total_size += vec->size();
    }

    //- Dummy resizing so we can pass *this to
    //- ListExpression::evaluate. Maybe do resizing externally?
    void resize_nocopy(const label)
    {}
    //- Assignment
    template<typename E>
    void operator=(const ListExpression<E, typename E::value_type>& expr)
    {
        expr.evaluate(*this);
    }

    value_type operator[](const label i) const
    {
        const auto ve(whichVecElem(i));
        return (*data[ve.first()])[ve.second()];
    }
    value_type& operator[](const label i)
    {    
        const auto ve(whichVecElem(i));
        auto& elems = *data[ve.first()];
        return elems[ve.second()];
    }
    label size() const { return total_size; }


    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef label difference_type;

    struct const_iterator
    {
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = this_type::difference_type;
        using value_type  = this_type::value_type;
        using pointer     = this_type::const_pointer;
        using reference   = this_type::const_reference;

        const this_type& base_; 
        label index_;

        const_iterator(const this_type& base, const label index)
        :
            base_(base),
            index_(index)
        {}
        const_iterator(const const_iterator& it)
        :
            base_(it.base_),
            index_(it.index_)
        {}
        const_iterator(const_iterator& it)
        :
            base_(it.base_),
            index_(it.index_)
        {}
        const_iterator() = delete;  //: base_(0) {}
        const_iterator& operator=(const_iterator& it) = default;

    // Indexing operators

        value_type operator*() const
        {
             return base_[index_];
        }
        //value_type operator->() const
        //{
        //    return *base_;
        //}
        value_type operator[](const difference_type i) const
        {
           return base_[index_+i];
        }

    // Increment/Decrement

        const_iterator& operator++()
        {
            ++index_;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator old(*this);
            ++index_;
            return old;
        }
        const_iterator& operator--()
        {
            --index_;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator old(*this);
            --index_;
            return old;
        }

    // Arithmetic operators

        const_iterator operator+(const difference_type i) const
        {
            return const_iterator(base_, index_+i);
        }
        const_iterator operator-(const difference_type i) const
        {
            return const_iterator(base_, index_-i);
        }
        difference_type operator-(const const_iterator& it) const
        {
            return (index_ - it.index_);
        }

    // Comparison operators

        bool operator==(const const_iterator& it) const
        {
            return (index_ == it.index_);
        }
        bool operator!=(const const_iterator& it) const
        {
            return (index_ != it.index_);
        }
        bool operator<(const const_iterator& it) const
        {
            return (index_ < it.index_);
        }
        bool operator<=(const const_iterator& it) const
        {
            return (index_ <= it.index_);
        }
        bool operator>(const const_iterator& it) const
        {
            return (index_ > it.index_);
        }
        bool operator>=(const const_iterator& it) const
        {
            return (index_ >= it.index_);
        }

    // Compound assignment 

        const_iterator& operator+=(const difference_type n)
        {
            index_ += n;
            return *this;
        }
        const_iterator& operator-=(const difference_type n)
        {
            index_ -= n;
            return *this;
        }
    };
    auto cbegin() const
    {
        return const_iterator(*this, 0);
    }
    auto cend() const
    {
        return const_iterator(*this, this->size());
    }

    struct iterator
    {
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = this_type::difference_type;
        using value_type  = this_type::value_type;
        using pointer     = this_type::pointer;
        using reference   = this_type::reference;

        this_type& base_; 
        label index_;

        iterator(this_type& base, const label index)
        :
            base_(base),
            index_(index)
        {}
        iterator(const iterator& it) = default;
        iterator() = delete;
        iterator& operator=(const iterator& it)
        {
            base_ = it.base_;
            index_ = it.index_;
            return *this;
        };

    // Indexing operators

        value_type& operator*()
        {
            return base_[index_];
        }
        value_type operator*() const
        {
            return base_[index_];
        }
        value_type operator->() const
        {
            return base_[index_];
        }
        value_type& operator->()
        {
            return base_[index_];
        }
        value_type operator[](const difference_type i) const
        {
           return base_[index_+i];
        }
        value_type& operator[](const difference_type i)
        {
           return base_[index_+i];
        }

    // Increment/Decrement

        iterator& operator++()
        {
            ++index_;
            return *this;
        }
        iterator operator++(int)
        {
            iterator old(*this);
            ++index_;
            return old;
        }
        iterator& operator--()
        {
            --index_;
            return *this;
        }
        iterator operator--(int)
        {
            iterator old(*this);
            --index_;
            return old;
        }

    // Arithmetic operators

        iterator operator+(const difference_type i) const
        {
            return iterator(base_, index_+i);
        }
        iterator operator-(const difference_type i) const
        {
            return iterator(base_, index_-i);
        }
        difference_type operator-(const iterator& it) const
        {
            return (index_ - it.index_);
        }

    // Comparison operators

        bool operator==(const iterator& it) const
        {
            return (index_ == it.index_);
        }
        bool operator!=(const iterator& it) const
        {
            return (index_ != it.index_);
        }
        bool operator<(const iterator& it) const
        {
            return (index_ < it.index_);
        }
        bool operator<=(const iterator& it) const
        {
            return (index_ <= it.index_);
        }
        bool operator>(const iterator& it) const
        {
            return (index_ > it.index_);
        }
        bool operator>=(const iterator& it) const
        {
            return (index_ >= it.index_);
        }

    // Compound assignment 

        iterator& operator+=(const difference_type n)
        {
            index_ += n;
            return *this;
        }
        iterator& operator-=(const difference_type n)
        {
            index_ -= n;
            return *this;
        }
    };
    auto begin()
    {
        auto& this_vals = *this;
        return iterator(this_vals, 0);
    }
    auto end()
    {
        auto& this_vals = *this;
        return iterator(this_vals, this->size());
    }
};


} // End namespace Expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
