from setuptools import setup

install_requires = [
    'six', 'setuptools >= 0.7', 'simplesam >= 0.0.3', 'tqdm', 'fastqp >= 0.2',
    'pandas', 'strandex', 'gffutils', 'pyfaidx', 'drmaa', 'twobitreader', 'requests'
]

setup(
    name='novartis-pisces',
    author='Matthew Shirley',
    author_email='matt_d.shirley@novartis.com',
    url='https://github.com/Novartis/pisces',
    description=
    'PISCES: pipeline for rapid transcript quantitation, genetic fingerprinting, and quality control assessment of RNAseq libraries',
    license='Apache',
    packages=['pisces'],
    include_package_data=True,
    install_requires=install_requires,
    use_scm_version={"local_scheme": "no-local-version"},
    setup_requires=['setuptools_scm'],
    entry_points={'console_scripts': ['pisces = pisces.cli:main', 'pisces_dependencies = pisces:install_dependencies']},
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Environment :: Console", "Intended Audience :: Science/Research",
        "Natural Language :: English", "Operating System :: Unix",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ])
