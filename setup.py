from setuptools import setup


def readme():
    with open("README.md", "r") as f:
        return f.read()

setup(name='irispreppy',
      version="0.9.4",
      url='https://github.com/OfAaron3/irispreppy',
      author='Aaron W. Peat',
      author_email='a.peat.1@research.gla.ac.uk',
      license='MIT',
      packages=['irispreppy.psf', 'irispreppy.radcal', 'irispreppy', 'irispreppy.radcal.responses'],
      install_requires=[
          'numpy',
          'tqdm',
          'astropy',
          'scipy',
          'beautifulsoup4'
      ],
      include_package_data=True,
      zip_safe=False,
      long_description=readme(),
      long_description_content_type="text/markdown",
)
