# based on https://realpython.com/pypi-publish-python-package
import pathlib
from setuptools import setup

from fast_fisher import __version__
from fast_fisher.fast_fisher_numba import cc

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / 'README.md').read_text()

if __name__ == "__main__":
    cc.compile()

    # This call to setup() does all the work
    setup(
        name='fast_fisher',
        version=__version__,
        description="Calculate Fisher's exact test very quickly.",
        long_description=README,
        long_description_content_type='text/markdown',
        url='https://github.com/mrtomrod/fast_fisher',
        author='Thomas Roder',
        author_email='roder.thomas@gmail.com',
        license='MIT',
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.10',
        ],
        packages=['fast_fisher'],
        install_requires=['numba'],
        ext_modules=[cc.distutils_extension()],
    )
