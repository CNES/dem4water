#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = []

test_requirements = []

setup(
    author="Benjamin Tardy",
    author_email="benjamin.tardy@csgroup.eu",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="dem4water estimate a law between altitude and surface to estimate water volume in dam.",
    entry_points={
        "console_scripts": [
            "dem4water=dem4water.cli:main",
        ],
    },
    install_requires=requirements,
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="dem4water",
    name="dem4water",
    packages=find_packages(include=["dem4water", "dem4water.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/tardyb/dem4water",
    version="0.1.0",
    zip_safe=False,
)
