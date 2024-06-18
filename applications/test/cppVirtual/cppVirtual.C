/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include <iostream>

/*---------------------------------------------------------------------------*\
                               Class a,b Declaration
\*---------------------------------------------------------------------------*/

class a
{
protected:
        int num_;
public:
        //- Constructor
        a(int num) : num_(num) {}

        //- Destructor
        ~a() {}

        //- Test a
        void test() const
        {
            std::cout << "a = " << num_ << std::endl;
        }
};

class b
:
    public a
{
public:
        //- Constructor
        b(int num) : a(num) {}

        //- Destructor
        ~b() {}

        //- Test b
        void test() const
        {
            std::cout << "b = " << num_ << std::endl;
        }
};

/*---------------------------------------------------------------------------*\
                              Class A,B Declaration
\*---------------------------------------------------------------------------*/

template<class ELEMENT>
class ptrListDummy
{
private:
        ELEMENT element_;
public:
        //- Null constructor
        ptrListDummy() : element_(NULL) {}

        //- Element constructor
        ptrListDummy(ELEMENT element) : element_(element) {}

        //- Set element
        bool set(ELEMENT element)
        {
            if (!element_)
            {
                element_ = element;
                return true;
            }
            else
            {
                std::cerr << "Element already exists!" << std::endl;
                return false;
            }
        }

        //- Get element
        ELEMENT& get()
        {
            return element_;
        }
};


class A
{
protected:
        int num_;
private:
        ptrListDummy<a*> dataA_;
public:
        //- Constructor
        A(int num, bool init = true)
        :
        num_(num),
        dataA_()
        {
            if (init) dataA_.set(new a(num));
        }

        virtual ~A() { if(dataA_.get()) delete dataA_.get(); };

        //- Test A
        void test() const
        {
            std::cout << "A = " << num_ << std::endl;
        }

        // Get data

            //- Static get data
            a& staticGetData()
            {
                return *dataA_.get();
            }

            //- Virtual get data
            virtual a& virtualGetData()
            {
                return *dataA_.get();
            }

        // Direct data test

            //- Static data test
            void staticDataTest()
            {
                return staticGetData().test();
            }

            //- Virtual data test
            virtual void virtualDataTest()
            {
                return virtualGetData().test();
            }

        // Direct base data test

            //- Static base data test
            void staticDataTestBase()
            {
                return staticGetData().test();
            }

            //- Static-virtual base data test
            void staticVirtualDataTestBase()
            {
                return virtualGetData().test();
            }

            //- Virtual base data test
            virtual void virtualDataTestBase()
            {
                return virtualGetData().test();
            }

            //- Virtual-static base data test
            virtual void virtualStaticDataTestBase()
            {
                return staticGetData().test();
            }
};


class B
:
    public A
{
private:
        ptrListDummy<b*> dataB_;
public:
        //- Constructor
        B(int num, bool init = true)
        :
        A(num, false),
        dataB_()
        {
            if (init) dataB_.set(new b(num));
        }

        virtual ~B() { if(dataB_.get()) delete dataB_.get(); };

        //- Test B
        void test() const
        {
            std::cout << "B = " << num_ << std::endl;
        }

        // Get data

            //- Static get data
            b& staticGetData()
            {
                return *dataB_.get();
            }

            //- Virtual get data
            virtual b& virtualGetData()
            {
                return *dataB_.get();
            }

        // Direct data test

            //- Static data test
            void staticDataTest()
            {
                return staticGetData().test();
            }

            //- Virtual data test
            virtual void virtualDataTest()
            {
                return virtualGetData().test();
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    std::cout << std::endl;
    std::cout << "################################# Basic test ###################################";
    std::cout << std::endl << std::endl;

        a my_a = a(0); std::cout << "a      => test()                          => ";
            my_a.test();
        b my_b = b(1); std::cout << "b      => test()                          => ";
            my_b.test();

    std::cout << std::endl;

        A my_A = A(2); std::cout << "A      => test()                          => ";
            my_A.test();
        B my_B = B(3); std::cout << "B      => test()                          => ";
            my_B.test();

    std::cout << std::endl;
    std::cout << "############################ Static get data test ##############################";
    std::cout << std::endl << std::endl;

        std::cout << "A      => staticGetData().test()          => ";
            my_A.staticGetData().test();
        std::cout << "B      => staticGetData().test()          => ";
            my_B.staticGetData().test();

    std::cout << std::endl;

        std::cout << "B      => staticGetData() => a => test()  => ";
            a& my_staticBa = my_B.staticGetData(); my_staticBa.test();
        std::cout << "B      => staticGetData() => b => test()  => ";
            b& my_staticBb = my_B.staticGetData(); my_staticBb.test();
        std::cout << "B => A => staticGetData().test()          => ";
//             A& my_staticBA = my_B; my_staticBA.staticGetData().test();
            std::cout << "SEGFAULT" << std::endl;

    std::cout << std::endl;
    std::cout << "########################### Static direct data test ############################";
    std::cout << std::endl << std::endl;

        std::cout << "A      => staticDataTest()                => ";
            my_A.staticDataTest();
        std::cout << "B      => staticDataTest()                => ";
            my_B.staticDataTest();

    std::cout << std::endl;

        std::cout << "A      => staticDataTestBase()            => ";
            my_A.staticDataTestBase();
        std::cout << "B      => staticDataTestBase()            => ";
//             my_B.staticDataTestBase();
            std::cout << "SEGFAULT" << std::endl;
        std::cout << "B      => staticVirtualDataTestBase()     => ";
            my_B.staticVirtualDataTestBase();

    std::cout << std::endl;
    std::cout << "############################ Virtual get data test #############################";
    std::cout << std::endl << std::endl;

        std::cout << "A      => virtualGetData().test()         => ";
            my_A.virtualGetData().test();
        std::cout << "B      => virtualGetData().test()         => ";
            my_B.virtualGetData().test();

    std::cout << std::endl;

        std::cout << "B      => virtualGetData() => a => test() => ";
            a& my_virtualBa = my_B.virtualGetData(); my_virtualBa.test();
        std::cout << "B      => virtualGetData() => b => test() => ";
            b& my_virtualBb = my_B.virtualGetData(); my_virtualBb.test();
        std::cout << "B => A => virtualGetData().test()         => ";
            A& my_virtualBA = my_B; my_virtualBA.virtualGetData().test();

    std::cout << std::endl;
    std::cout << "########################### Virtual direct data test ###########################";
    std::cout << std::endl << std::endl;

        std::cout << "A      => virtualDataTest()               => ";
            my_A.virtualDataTest();
        std::cout << "B      => virtualDataTest()               => ";
            my_B.virtualDataTest();

    std::cout << std::endl;

        std::cout << "A      => virtualDataTestBase()           => ";
            my_A.virtualDataTestBase();
        std::cout << "B      => virtualDataTestBase()           => ";
            my_B.virtualDataTestBase();
        std::cout << "B      => virtualStaticDataTestBase()     => ";
//             my_B.virtualStaticDataTestBase();
            std::cout << "SEGFAULT" << std::endl;

    std::cout << std::endl;

    std::cout << "################################################################################";
    std::cout << std::endl << std::endl;

    return 0;
}

// ************************************************************************* //

// ################################# Basic test ###################################
//
// a      => test()                          => a = 0
// b      => test()                          => b = 1
//
// A      => test()                          => A = 2
// B      => test()                          => B = 3
//
// ############################ Static get data test ##############################
//
// A      => staticGetData().test()          => a = 2
// B      => staticGetData().test()          => b = 3
//
// B      => staticGetData() => a => test()  => a = 3
// B      => staticGetData() => b => test()  => b = 3
// B => A => staticGetData().test()          => SEGFAULT
//
// ########################### Static direct data test ############################
//
// A      => staticDataTest()                => a = 2
// B      => staticDataTest()                => b = 3
//
// A      => staticDataTestBase()            => a = 2
// B      => staticDataTestBase()            => SEGFAULT
// B      => staticVirtualDataTestBase()     => a = 3
//
// ############################ Virtual get data test #############################
//
// A      => virtualGetData().test()         => a = 2
// B      => virtualGetData().test()         => b = 3
//
// B      => virtualGetData() => a => test() => a = 3
// B      => virtualGetData() => b => test() => b = 3
// B => A => virtualGetData().test()         => a = 3
//
// ########################### Virtual direct data test ###########################
//
// A      => virtualDataTest()               => a = 2
// B      => virtualDataTest()               => b = 3
//
// A      => virtualDataTestBase()           => a = 2
// B      => virtualDataTestBase()           => a = 3
// B      => virtualStaticDataTestBase()     => SEGFAULT
//
// ################################################################################

// ************************************************************************* //
