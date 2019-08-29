#!/usr/bin/env python3
from setuptools import setup, find_packages
from glob import glob
import versioneer


desc = """
pysolarcalc: Simulate temp/humidiy/light from arbirtary points on the earth for
environment-mimicking growth cabinets
"""

with open("requirements.txt") as fh:
    install_requires = [req.strip() for req in fh]

test_requires = [
    "pytest",
]


setup(
    name="pysolarcalc",
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=install_requires,
    tests_require=test_requires,
    description=desc,
    entry_points='''
        [console_scripts]
        pysolarcalc=solarcalc.cli:main
    ''',
    author="Kevin Murray",
    author_email="foss@kdmurray.id.au",
    url="https://github.com/appf-anu/pysolarcalc",
    keywords=["science", "solarcalc"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    test_suite="tests",
    )

