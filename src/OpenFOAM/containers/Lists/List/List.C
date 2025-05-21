/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "FixedList.H"
#include "PtrList.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::List<T>::resize_copy(const label count, const label len)
{
    // Only a limited number of internal size checks.
    // Caller knows what they are doing.

    if (FOAM_LIKELY(len > 0))
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        const label oldLen = this->size_;
        const label overlap = Foam::min(count, len);
        // Extra safety, not currently necessary:
        // const label overlap = Foam::min(Foam::min(count, oldLen), len);

        T* old = this->v_;

        if (overlap > 0)
        {
            // Recover overlapping content when resizing

            this->size_ = len;
            this->v_ = ListPolicy::allocate<T>(len);

            // Can dispatch with
            // - std::execution::par_unseq
            // - std::execution::unseq
            std::move(old, (old + overlap), this->v_);

            ListPolicy::deallocate(old, oldLen);
        }
        else
        {
            // No overlapping content
            ListPolicy::deallocate(old, oldLen);

            this->size_ = len;
            this->v_ = ListPolicy::allocate<T>(len);
        }
    }
    else
    {
        // Or only #ifdef FULLDEBUG
        if (FOAM_UNLIKELY(len < 0))
        {
            FatalErrorInFunction
                << "bad size " << len
                << abort(FatalError);
        }
        // #endif

        clear();
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(const label len)
:
    UList<T>(nullptr, len)
{
    if (FOAM_UNLIKELY(len < 0))
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len > 0)
    {
        doAlloc();
    }
}


template<class T>
Foam::List<T>::List(const label len, const T& val)
:
    UList<T>(nullptr, len)
{
    if (FOAM_UNLIKELY(len < 0))
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len > 0)
    {
        doAlloc();
        UList<T>::operator=(val);
    }
}


template<class T>
Foam::List<T>::List(const label len, Foam::zero)
:
    UList<T>(nullptr, len)
{
    if (FOAM_UNLIKELY(len < 0))
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len > 0)
    {
        doAlloc();
        UList<T>::operator=(Foam::zero{});
    }
}


template<class T>
Foam::List<T>::List(Foam::one, const T& val)
:
    UList<T>(ListPolicy::allocate<T>(1), 1)
{
    this->v_[0] = val;
}


template<class T>
Foam::List<T>::List(Foam::one, T&& val)
:
    UList<T>(ListPolicy::allocate<T>(1), 1)
{
    this->v_[0] = std::move(val);
}


template<class T>
Foam::List<T>::List(Foam::one, Foam::zero)
:
    UList<T>(ListPolicy::allocate<T>(1), 1)
{
    this->v_[0] = Foam::zero{};
}


template<class T>
Foam::List<T>::List(const UList<T>& list)
:
    UList<T>(nullptr, list.size_)
{
    if (this->size_ > 0)
    {
        doAlloc();
        UList<T>::deepCopy(list);
    }
}


template<class T>
Foam::List<T>::List(const List<T>& list)
:
    UList<T>(nullptr, list.size_)
{
    if (this->size_ > 0)
    {
        doAlloc();
        UList<T>::deepCopy(list);
    }
}


template<class T>
Foam::List<T>::List(List<T>& list, bool reuse)
:
    UList<T>(nullptr, list.size_)
{
    if (reuse)
    {
        // Steal content
        this->v_ = list.v_;
        list.v_ = nullptr;
        list.size_ = 0;
    }
    else if (this->size_ > 0)
    {
        doAlloc();
        UList<T>::deepCopy(list);
    }
}


template<class T>
Foam::List<T>::List(const UList<T>& list, const labelUList& indices)
:
    UList<T>(nullptr, indices.size())
{
    doAlloc();
    copyList(list, indices);  // <- deepCopy()
}


template<class T>
template<unsigned N>
Foam::List<T>::List
(
    const UList<T>& list,
    const FixedList<label,N>& indices
)
:
    UList<T>(nullptr, indices.size())
{
    doAlloc();
    copyList(list, indices);  // <- deepCopy()
}


template<class T>
template<unsigned N>
Foam::List<T>::List(const FixedList<T, N>& list)
:
    List<T>(list.begin(), list.end(), list.size())
{}


template<class T>
Foam::List<T>::List(const PtrList<T>& list)
:
    UList<T>(nullptr, list.size())
{
    doAlloc();
    copyList(list);
}


template<class T>
template<class Addr>
Foam::List<T>::List(const IndirectListBase<T, Addr>& list)
:
    UList<T>(nullptr, list.size())
{
    if (this->size_ > 0)
    {
        doAlloc();
        UList<T>::deepCopy(list);
    }
}


template<class T>
Foam::List<T>::List(std::initializer_list<T> list)
:
    List<T>(list.begin(), list.end(), list.size())
{}


template<class T>
Foam::List<T>::List(List<T>&& list) noexcept
:
    UList<T>(list.data(), list.size())
{
    list.size_ = 0;
    list.v_ = nullptr;
}


template<class T>
template<int SizeMin>
Foam::List<T>::List(DynamicList<T, SizeMin>&& list)
:
    UList<T>()
{
    transfer(list);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::~List()
{
    //TODO? May need to verify that size is accurate (for correct alignment)
    ListPolicy::deallocate(this->v_, this->size_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::resize(const label len, const T& val)
{
    const label oldLen = this->size_;

    if (oldLen == len)
    {
        return;
    }

    this->resize_copy(oldLen, len);

    // Fill trailing part with new values
    if (oldLen < this->size_)
    {
        std::fill
        (
            (this->v_ + oldLen), (this->v_ + this->size_), val
        );
    }
}


template<class T>
void Foam::List<T>::transfer(List<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // Clear and swap
    clear();
    this->size_ = list.size_;
    this->v_ = list.v_;

    list.size_ = 0;
    list.v_ = nullptr;
}


template<class T>
template<int SizeMin>
void Foam::List<T>::transfer(DynamicList<T, SizeMin>& list)
{
    // Remove existing contents before anything else.
    clear();

    // Shrink the allocated space to the number of elements used
    list.shrink_to_fit();
    transfer(static_cast<List<T>&>(list));
    list.clearStorage();  // Deletion, capacity=0 etc.
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::operator=(const UList<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    reAlloc(list.size_);

    if (this->size_ > 0)
    {
        UList<T>::deepCopy(list);
    }
}


template<class T>
void Foam::List<T>::operator=(const List<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    reAlloc(list.size_);

    if (this->size_ > 0)
    {
        UList<T>::deepCopy(list);
    }
}


template<class T>
template<unsigned N>
void Foam::List<T>::operator=(const FixedList<T, N>& list)
{
    reAlloc(list.size());

    std::copy(list.begin(), list.end(), this->v_);
}


template<class T>
template<class Addr>
void Foam::List<T>::operator=(const IndirectListBase<T, Addr>& list)
{
    reAlloc(list.size());
    UList<T>::deepCopy(list);
}


template<class T>
void Foam::List<T>::operator=(std::initializer_list<T> list)
{
    reAlloc(list.size());

    std::copy(list.begin(), list.end(), this->v_);
}


template<class T>
void Foam::List<T>::operator=(List<T>&& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    transfer(list);
}


template<class T>
template<int SizeMin>
void Foam::List<T>::operator=(DynamicList<T, SizeMin>&& list)
{
    transfer(list);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
Foam::labelList Foam::sortedOrder
(
    const UList<T>& list
)
{
    labelList order;
    Foam::sortedOrder(list, order, typename UList<T>::less(list));
    return order;
}


template<class T>
void Foam::sortedOrder
(
    const UList<T>& list,
    labelList& order
)
{
    Foam::sortedOrder(list, order, typename UList<T>::less(list));
}


template<class T, class ListComparePredicate>
void Foam::sortedOrder
(
    const UList<T>& list,
    labelList& order,
    const ListComparePredicate& comp
)
{
    // List lengths must be identical. Old content is overwritten
    order.resize_nocopy(list.size());

    Foam::identity(order, 0);
    Foam::stableSort(order, comp);
}


// ************************************************************************* //
