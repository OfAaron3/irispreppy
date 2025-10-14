from setuptools import setup


def readme():
    with open("README.md", "r") as f:
        return f.read()

setup(name='irispreppy',
      description="For radiometric calibration and point spread function deconvolution of IRIS spsetcrograph data",
      setup_requires=["setuptools_scm"],
      use_scm_version=True,
      url='https://github.com/OfAaron3/irispreppy',
      author='Aaron W. Peat',
      author_email='aaron.peat@uwr.edu.pl',
      license='MIT',
      packages=['irispreppy.psf', 'irispreppy.radcal', 'irispreppy', 'irispreppy.radcal.responses'],
      install_requires=[
          'numpy',
          'astropy',
          'scipy',
          'beautifulsoup4',
          'weno4'
      ],
      include_package_data=True,
      zip_safe=False,
      long_description=readme(),
      long_description_content_type="text/markdown",
)
