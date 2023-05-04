from setuptools import setup

setup(name='SSAxgeo',
      version='1.0.dev',
      description='''
      Protein secondary structure elements assignment based on differential geometry
      ''',
      url='https://github.com/labstructbioinf/SSAxgeo',
      author='"Antonio Marinho',
      author_email='amarinhosn@pm.me',
      packages=['ssaxgeo'],
      classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
          ],
      scripts=['scripts/ssaxgeo', 'scripts/getSampleOfClstrPDB', 'scripts/computePDBxgeo'],
      install_requires=['hdbscan', 'pandas', 'numpy', 
                        'matplotlib', 'seaborn', 'scipy', 'localpdb'],
      zip_safe=False)