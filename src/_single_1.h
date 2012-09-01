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

	PyReadGraphLoc &operator=(PyReadGraphLoc const &x) {
		if (this != &x) {
			graph_id = x.graph_id;
			locs = x.locs;
		}
		return *this;
	}

	PyReadGraphLoc(PyReadGraphLoc const &x) :
			graph_id(x.graph_id), locs(x.locs) {
	}

	~PyReadGraphLoc() {
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

}

#endif // _SINGLE_1_H_
