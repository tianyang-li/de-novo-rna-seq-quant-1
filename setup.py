#  Copyright (C) 2012 Tianyang Li
#  tmy1018@gmail.com
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License


from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


cflags = ['-Wall', '-Wextra', '-Wconversion',
          #'-ansi', '-pedantic',
          '-Wall',
          '-W', '-Wmissing-prototypes',
          '-Wstrict-prototypes',
          # '-Wpointer-arith', '-Wcast-qual', '-Wcast-align',
          '-Wwrite-strings', '-Wnested-externs',
          '-fshort-enums',
          '-fno-common',
          #'-Wshadow',
          ]


setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("single_1",
                           ["src/single_1.pyx"],
                           libraries=[],
                           extra_link_args=['-fPIC'] + cflags,
                           extra_compile_args=['-fPIC'] + cflags,
                           extra_objects=["lib/quant.a", "lib/libgsl.a"],
                           language="c++",
                           include_dirs=["src", "gsl-1.15", "boost_1_51_0"],
                           library_dirs=["lib"]
                           )
                 ]
      )




