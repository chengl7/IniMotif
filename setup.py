from distutils.core import setup

setup(
      name = "IniMotif",
      version = "0.1.0",
      python_requires='>=3.0',
      install_requires = [
                          'numpy',
                          'scipy',
                          'matplotlib',
                          'biopython',
                          'yattag',
                          'adjustText',
                          'seaborn',
                          'requests',
                          ],
      classifiers = [
                     "Programming Language :: Python :: 3",
                     "Operating System :: OS Independent",
                     ],
      )
