from setuptools import setup

setup(name='SSAxgeo',
      version='1.0.dev',
      description='''
      Workflows for SARS-CoV-2 genome Assembly at FioCruz/IAM
      ''',
      url='https://github.com/dezordi/ViralFlow/',
      author='"Antonio Marinho',
      author_email='amarinhosn@pm.me',
      packages=['ssaxgeo'],
      classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
          ],
      scripts=['scripts/ssaxgeo'],
      install_requires=['hdbscan', 'pandas', 'numpy', 
                        'matplotlib', 'seaborn', 'scipy'],
      zip_safe=False)