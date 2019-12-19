from distutils.core import setup

setup(
      name = "IniMotif",
      version = "0.1.0",
      packages = setuptools.find_packages(),
      scripts = [
                 'inimotif_gui.py',
                 'dna_logo.py',
                 'inimotif_core.py',
                 'inimotif_main.py',
                 'inimotif_util.py',
                 'windows.py',
                 'chipWinExtract.py',
                 'enaFastqFetch.py',
                 'xmlparser.py',
                 'reportwriter.py',
                 ],
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
      data_files=[('GUIgraphics', ['GUIgraphics/logo.png', 'GUIgraphics/masker.png'])],
      include_package_data=True,
      )
