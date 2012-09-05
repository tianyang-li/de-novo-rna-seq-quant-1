/*
 * Tianyang Li 2012
 * tmy1018@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MCMC_0_H_
#define _MCMC_0_H_

#include "_graph_seq_0.h"

namespace _mcmc_0 {

using _graph_seq_0::SpliceGraph;
using _graph_seq_0::GraphReads;

// calculate the probability that a read
// is from a the transcripts (isoforms)
template<class RNodeLoc>
ldbl ReadFromTransProb();

}

#endif // _MCMC_0_H_
