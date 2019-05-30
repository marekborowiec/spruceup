from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as fr:
    requirements = fr.read()

setup(
    name='example-pkg-your-username',
    version='2019.1.0',
    author='Marek Borowiec',
    author_email='petiolus@gmail.com',
    description='A module for lexible identification, visualization, and removal of outliers from large multiple sequence alignments',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/marekborowiec/seq-spruceup',
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research'
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        'Environment :: Console',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
)
