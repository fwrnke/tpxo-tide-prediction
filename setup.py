import os
import setuptools

def read(fname, mode='r', encoding='utf-8'):
    """Wrap open() function."""
    return open(os.path.join(os.path.dirname(__file__), 'tpxo_tide_prediction', fname), mode=mode, encoding=encoding)

with open('README.md') as f:
    long_description = f.read()

VERSION = {}
with read('_version.py') as fp:
    exec(fp.read(), VERSION)

setuptools.setup(
    name='tpxo-tide-prediction',
    version=VERSION['__version__'],
    author='Fynn Warnke',
    author_email='fwrnke@mailbox.org',
    description='Predict tidal elevation based on TPXO9-atlas model',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GNU GPLv3+',
    url='https://github.com/fwrnke/tpxo-tide-prediction',
    project_urls={
        'Bug Tracker': 'https://github.com/fwrnke/tpxo-tide-prediction/issues',
    },
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering'
    ],
    packages=setuptools.find_packages(exclude=('tests')),
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.16.5',
        'xarray>=0.20.0',
        'scipy>=1.7.0',
    ],
    tests_require=['pytest', 'pytest-allclose'],
    entry_points={
    'console_scripts': [
        'predict_tide=tpxo_tide_prediction.predict_tide:main',
        ],
    },
)
