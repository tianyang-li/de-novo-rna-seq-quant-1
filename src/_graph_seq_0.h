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

/*
 * some general stuff
 * for working with sequences, graphs, etc.
 */

#ifndef _GRAPH_SEQ_0_H_
#define _GRAPH_SEQ_0_H_

#include "_misc_0.h"

namespace _graph_seq_0 {

/*
 * assumes that we already know which graph it is
 */
class SeqLoc {
public:
	SeqLoc(int node_id_, uint start_, uint end_) :
			node_id(node_id_), start(start_), end(end_) {
	}

	SeqLoc &operator=(SeqLoc const &x) {
		if (this != &x) {
			node_id = x.node_id;
			start = x.start;
			end = x.end;
		}
		return *this;
	}

	SeqLoc(SeqLoc const &x) :
			node_id(x.node_id), start(x.start), end(x.end) {
	}

	~SeqLoc() {
	}

	int node_id;
	uint start;
	uint end;
};

}

#endif // _GRAPH_SEQ_0_H_
