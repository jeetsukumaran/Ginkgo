///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical and Evolution Simulator.
//
// Copyright 2009 Jeet Sukumaran and Mark T. Holder.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#if !defined(GINGKO_RANDOM_H)
#define GINGKO_RANDOM_H

namespace gingko {

/**
 * Encapsulates generation of random variates drawn from various distributions.
 */
class RandomNumberGenerator {

    public:
        RandomNumberGenerator();
        RandomNumberGenerator(unsigned long seed);
        void set_seed(unsigned long seed);
 
        float uniform_real();                   //! Uniform random real variate in [0, 1).
        long  uniform_int(int a, int b);        //! Uniform random integer variate in [a, b].
        float standard_normal();                //! Gaussian random variate with mean of 0 and std. dev. of 1.
        float normal(float mean, float sd);     //! Gaussian random variate with given mean and std. dev.        
        unsigned int poisson(float rate);       //! Poisson random variate with given hazard parameter.
        
        //! Returns an element selected with uniform probability from a
        //! collection.
        template <typename T>
        inline typename T::value_type& select(T& collection) {
            return collection[this->uniform_int(0, collection.size()-1)];
        }
        
        //! Returns one of two values passed as arguments with uniform
        //! random probability.
        template <typename T>        
        inline T& select(T& a, T& b) {
            if (this->uniform_real() < 0.5) {
                return a;
            } else {
                return b;
            }
        }        
        
    private:
        unsigned long seed_;                    //! seed for the underlying rng
};


} // gingko namespace

#endif