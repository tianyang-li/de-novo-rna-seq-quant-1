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

/*
 * assumes that the read is already known
 */
class PyReadGraphLoc {
public:
	PyReadGraphLoc(uint graph_id_) :
			graph_id(graph_id_) {
	}

	uint graph_id;
	std::vector<_graph_seq_0::PyReadNodeLoc> locs;
};

class PyReadInGraph {
public:
	PyReadInGraph() {
	}

	uint read_id;
	std::vector<PyReadGraphLoc> graph_locs;
};

/*
 * used to interface with python
 */
void _get_isoforms(std::vector<_graph_seq_0::PyGraph> *py_graphs,
		std::vector<PyReadInGraph> *py_reads,
		std::vector<_graph_seq_0::Isoform> *isoforms);

/*
 * keep the reads aligned to a particular graph
 */

/*
 * given read's @read_id
 * the alignment is at
 * py_reads[read_id].graph_locs[graph_index].locs[align_index]
 */
class ReadIndex {
public:
	ReadIndex(uint read_id_, uint graph_index_, uint align_index_) :
			read_id(read_id_), graph_index(graph_index_), align_index(
					align_index_) {
	}

	uint read_id;
	uint graph_index;
	uint align_index;
};

class GraphReads {
public:
	uint graph_id;
	std::vector<ReadIndex> reads;
};

}

#endif // _SINGLE_1_H_
