from setuptools import setup

setup(name='irispreppy',
      version=0.9,
      url='https://github.com/OfAaron3/irispreppy',
      author='Aaron W. Peat',
      author_email='a.peat.1@research.gla.ac.uk',
      license='MIT',
      packages=['irispreppy.psf', 'irispreppy.radcal'],
      install_requires=[
          'numpy',
          'tqdm',
          'astropy',
          'scipy',
          'bs4'
      ],
      include_package_data=True,
      zip_safe=False
)
