from setuptools import setup

setup(name='irispreppy',
      version=0.9,
      url='https://github.com/OfAaron3/irispreppy',
      author='Aaron W. Peat',
      author_email='a.peat.1@research.gla.ac.uk',
      license='MIT',
      py_modules=['deconvolve', 'radcal'],
      install_requires=[
          'numpy',
          'tqdm',
          'astropy',
          'scipy',
          'bs4'
      ],
      zip_safe=False
)