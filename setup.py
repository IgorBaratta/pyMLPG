from setuptools import setup

setup(name='pymlpg',
      version='0.0.1',
      description='Python package for the Meshless Local Petrov Galerkin method',
      url='https://github.com/IgorBaratta/pyMLPG',
      author='Igor Baratta',
      author_email='igorbaratta@gmail.com',
      license='MIT',
      packages=['pymlpg'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib'
      ],
      zip_safe=False)