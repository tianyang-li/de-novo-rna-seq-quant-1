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

// commonly used exceptions
#ifndef _ERR_0_H_
#define _ERR_0_H_

#include <exception>

namespace _err_0 {

using std::exception;

#ifdef DEBUG
class IteratorEndError: public exception {
public:
	inline virtual char const * what() const throw () {
		return "IteratorEndError:\n iterator didn't reach end as expected!\n";
	}
};
#endif

}

#endif // _ERR_0_H_
