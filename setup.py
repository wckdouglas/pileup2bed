from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy as np

setup(
    name='pileup2bed',
    version='0.1',
    description='parsing mpileup result to bed file wiwht base count',
    url='',
    author='Douglas Wu',
    author_email='wckdouglas@gmail.com',
    license='MIT',
    packages=['pileup2bed'],
    zip_safe=False,
    scripts = ['bin/pileup_to_bed.py'],
    ext_modules = cythonize([Extension('pileup2bed.parsing_pileup',
                                       ['pileup2bed/parsing_pileup.pyx'],
                                       include_dirs = [np.get_include()])]),
    cmdclass = {'build_ext': build_ext}
)
