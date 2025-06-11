/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseSystem.H"
#include "relativePhaseInterface.H"
//#include "relativeDisplacedPhaseInterface.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class wordListAndType
{
    public:

        wordList wl;
        Type t;

        wordListAndType()
        {}

        wordListAndType(Istream& is)
        :
            wl(is),
            t(is)
        {}
};

template<class Type>
inline Istream& operator>>(Istream& is, wordListAndType<Type>& wlat)
{
    return is >> wlat.wl >> wlat.t;
}

template<class Type>
inline Ostream& operator<<(Ostream& os, const wordListAndType<Type>& wlat)
{
    return os << wlat.wl << wlat.t;
}

template<class Type>
inline bool operator==
(
    const wordListAndType<Type>& a,
    const wordListAndType<Type>& b
)
{
    return a.wl == b.wl && a.t == b.t;
}

template<class Type>
inline bool operator!=
(
    const wordListAndType<Type>& a,
    const wordListAndType<Type>& b
)
{
    return !(a == b);
}

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class ModelType>
void Foam::phaseSystem::generateInterfacialModels
(
    const dictionary& dict,
    HashTable
    <
        autoPtr<ModelType>,
        phaseInterfaceKey,
        phaseInterfaceKey::hash
    >& models
) const
{
    // Construct lists of interfaces and models
    PtrList<phaseInterface> listInterfaces;
    PtrList<ModelType> listModels;

    // Set size of PtrLists. Resize could be necessary at the end.
    listModels.setSize(dict.size());
    listInterfaces.setSize(dict.size());

    label ptrIter = 0;

    forAllConstIter(dictionary, dict, iter)
    {
        // Get the model sub dictionary and its associated interface
        const dictionary& modelDict = iter().dict();
        autoPtr<phaseInterface> modelInterfacePtr =
            phaseInterface::New(*this, iter().keyword());

        // Cast the interface down to the first specified type possible
        autoPtr<phaseInterface> interfacePtr;
        
        if (isA<relativePhaseInterface>(modelInterfacePtr()))
        {
            interfacePtr.set
            (
                new relativePhaseInterface
                (
                    refCast<relativePhaseInterface>(modelInterfacePtr())
                )
            );
        }
        else if (isA<phaseInterface>(modelInterfacePtr()))
        {
            interfacePtr.set
            (
                new phaseInterface
                (
                    refCast<phaseInterface>(modelInterfacePtr())
                )
            );
        }

        if (!interfacePtr.valid())
        {
            FatalErrorInFunction
                << "Interface " << modelInterfacePtr->name()
                << " is not of suitable type for construction of a "
                << ModelType::typeName
                << exit(FatalError);
        }

        // If constructing for a specific interface then combine with this
        // interface. This ensures interface information propagates through
        // hierarchical model generation.
        //if (notNull(interface))
        //{
        //    interfacePtr = phaseInterface::New(interface, interfacePtr());
        //}

        listInterfaces.set(ptrIter, interfacePtr.ptr());
        listModels.set(ptrIter, ModelType::New(modelDict, listInterfaces[ptrIter]));
        ptrIter++;
    }

    // Transfer to a keyed table
    forAll(listInterfaces, i)
    {
        models.insert(listInterfaces[i], listModels.set(i, nullptr));
    }
}


// ************************************************************************* //
