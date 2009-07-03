///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical Evolution Simulator.
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
 * Encapsulates generation of random number from various distribution and 
 * various ranges.
 */
class RandomNumberGenerator {

    public:
    
        /** Constructs a RNG seeded with current time. */
        RandomNumberGenerator();
        
        /** Constructs a RNG seeded with current time. */
        
        /** Constructs a RNG seeded with given seed. */
        RandomNumberGenerator(unsigned long seed);
        
        /** Returns current seed. */
        unsigned long get_seed() const;        
        
        /** Explicitly sets the RNG seed. */
        void set_seed(unsigned long seed);                
 
        /**
         * Returns a uniform random real variate in [0,1).
         * @return   uniform random real variate in [0,1)
         */
        float uniform_01();
        
        /**
         * Returns a uniform random integer in [a, b].
         * @param   lower-bound of range
         * @param   upper-bound of range
         * @return  uniform random integer in [a, b]
         */        
        long  uniform_int(int a, int b);
        
        /**
         * Returns Gaussian random variate with mean of 0 and standard
         * deviation of 1.
         * @return   random variate with mean of 0 and standard deviation of 1
         */          
        float standard_normal(); 
        
        /**
         * Returns Gaussian random variate with mean of <code>mean</code> and 
         * standard deviation of <code>sd</code>.
         * @return   random variate with mean of <code>mean</code> and standard 
         *           deviation of <code>sd</code>
         */        
        float normal(float mean, float sd);
        
        /**
         * Returns a Poisson-distributed random variate with given rate.
         * @param   rate    rate for the Poisson process
         * @return          random variate from a Poisson distribution with given
         *                  rate
         */
        unsigned int poisson(float rate);
        
        /**
         * Returns an element selected with uniform random probability from
         * given universe of elements.
         * @param   collection  universe of elements from which to sample
         * @return              random element from collection
         */
        template <typename T>
        inline typename T::value_type& select(T& collection) {
            return collection[this->uniform_int(0, collection.size()-1)];
        }
        
        /**
         * Returns one of two arguments passed to it.
         * @param   a   the first candidate to be returned
         * @param   b   the second candidate to be returned
         * @return      either <code>a</code> or <code>b</code>, selected at
         *              random
         */
        template <typename T>        
        inline T& select(T& a, T& b) {
            if (this->uniform_01() < 0.5) {
                return a;
            } else {
                return b;
            }
        }        
        
    private:
        /** random number generator seed */
        unsigned long seed_;                    //! seed for the underlying rng
};

/**
 * Random function used as function object by STL algorithm.
 */
class RandomPointer {     
    public:     
        RandomPointer(RandomNumberGenerator& rng);            
        std::ptrdiff_t operator() (std::ptrdiff_t max);
    private:
        RandomNumberGenerator&  rng_;
};

} // gingko namespace

#endif
