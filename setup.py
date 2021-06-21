"""Setup for the locCSN package."""
with open('README.md') as f:
    README = f.read()

import setuptools

setuptools.setup(
    author="Xuran Wang",
    author_email="xuranw@andrew.cmu.edu",
    name='locCSN',
    license="MIT",
    description='locCSN is a python package for local cell specific networks.',
    version='0.0.11',
    url='https://github.com/xuranw/locCSN',
    packages=setuptools.find_packages(),
    python_requires=">=3.5",
    install_requires=['numpy', 'pandas', 'scipy', 'matplotlib', 'scanpy', 'joblib'],
    classifiers=[
        # Trove classifiers
        # (https://pypi.python.org/pypi?%3Aaction=list_classifiers)
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Intended Audience :: Developers',
    ],
)