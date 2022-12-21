# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='CEMCodesHub',
    version='0.1.0',
    description='A GUI for accelerating cavity design.',
    long_description=readme,
    author='Sosoho-Abasi Udongwo',
    author_email='numurho@gmail.com',
    url=r'https://github.com/Dark-Elektron/CEMCodesHub',
    license=license,
    python_requires='>=3.0, <4',
    install_requires=[
        'pandas==0.23.3',
        'numpy>=1.14.5',
        'matplotlib>=2.2.0'
    ],
    packages=find_packages(exclude=('tests', 'docs')),
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Topic :: Scientific/Engineering :: Physics",
        ],
)
