from setuptools import setup

setup(name='irispreppy',
      version="0.9.1",
      url='https://github.com/OfAaron3/irispreppy',
      author='Aaron W. Peat',
      author_email='a.peat.1@research.gla.ac.uk',
      license='MIT',
      packages=['irispreppy.psf', 'irispreppy.radcal', 'irispreppy', 'irispreppy.radcal.responses'],
      install_requires=[
          'numpy',
          'tqdm',
          'astropy<=4.2.0',
          'scipy',
          'beatifulsoup4'
      ],
      include_package_data=True,
      zip_safe=False
)
