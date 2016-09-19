from distutils.core import setup, Extension
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError("Requires cython to "
            "be installed before running setup.py (pip install cython)")
try:
    import numpy as np
except ImportError:
    raise ImportError("Requires numpy to "
            "be installed before running setup.py (pip install numpy)")

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
    install_requires=[
          'cython',
          'numpy'
      ],
    cmdclass = {'build_ext': build_ext}
)
