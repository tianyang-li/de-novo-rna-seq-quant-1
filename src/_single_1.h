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

#ifndef _SINGLE_1_H_
#define _SINGLE_1_H_

#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "_misc_0.h"
#include "_graph_seq_0.h"

namespace _single_1 {

using _graph_seq_0::ReadInGraph;
using _graph_seq_0::ReadIndex;

/*
 * assumes that the read is already known
 */
class SingleGraphLoc {
public:
	SingleGraphLoc(uint graph_id_) :
			graph_id(graph_id_) {
	}

	uint graph_id;
	std::vector<_graph_seq_0::ReadNodeLoc> locs;
};

/*
 * used to interface with python
 */
void _get_isoforms(std::vector<_graph_seq_0::PyGraph> *py_graphs,
		std::vector<ReadInGraph<SingleGraphLoc> > *py_reads,
		std::vector<_graph_seq_0::Fasta> *isoforms);

/*
 * keep the reads aligned to a particular graph
 */
class GraphReads {
public:
	uint graph_id;
	std::vector<ReadIndex> reads;
};

}

#endif // _SINGLE_1_H_
