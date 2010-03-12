///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
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

#if !defined(GINKGO_RANDOM_H)
#define GINKGO_RANDOM_H

#include <vector>

namespace ginkgo {

///////////////////////////////////////////////////////////////////////////////
// RandomNumberGenerator

/**
 * Encapsulates generation of random number from various distribution and
 * various ranges.
 */
class RandomNumberGenerator {

    public:

        // default destructor
        ~RandomNumberGenerator() {}

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
         * Returns an index value with probability equal to a vector of weights
         * passed as an argument. For example, given {2.0, 2.0, 3.0, 3.0}, this
         * will return:
         *       0 with probability 0.2
         *       1 with probability 0.2
         *       2 with probability 0.3
         *       3 with probability 0.3
         */
        unsigned int weighted_index_choice(const std::vector<float>& weights);

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
         * Returns a pointer to element selected with uniform random probability from
         * given universe of elements.
         * @param   collection  universe of elements from which to sample
         * @return              pointer to random element from collection
         */
        template <typename T>
        inline typename T::value_type* select_ptr(T& collection) {
            return &collection[this->uniform_int(0, collection.size()-1)];
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


    ///////////////////////////////////////////////////////////////////////////
    // Singleton infrastructure

    public:
        static RandomNumberGenerator& get_instance() {
            return RandomNumberGenerator::instance_;
        }

    private:
        static RandomNumberGenerator instance_;
        RandomNumberGenerator() {}
        RandomNumberGenerator(const RandomNumberGenerator &);
        RandomNumberGenerator & operator=(const RandomNumberGenerator &);

};

///////////////////////////////////////////////////////////////////////////////
// RandomPointer
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


} // ginkgo namespace



#endif
