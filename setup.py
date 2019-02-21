from setuptools import setup
from distutils.core import setup, Extension

gssw_extension = Extension('_gssw', sources=['pygssw/gssw_wrap.c', 'pygssw/gssw.c'], )


setup(name='pygssw',
      version='0.0.1',
      description='Wrapper around GSSW, making it possible to align reads to graphs from python using GSSW',
      url='',
      author='Ivar Grytten',
      author_email='ivargry@ifi.uio.no',
      license='MIT',
      zip_safe=False,
      install_requires=[],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      ext_modules=[gssw_extension],
      package_data={'': ['pygssw/gssw.h']}

)