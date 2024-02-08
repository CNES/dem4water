#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open("README.md", encoding="utf-8") as readme_file:
    readme = readme_file.read()

requirements = [
    "geopandas",
    "scipy",
    "numpy",
    "matplotlib",
    "flake8",
    "gdal",
    "rasterio",
    "bmi_topography",
]


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
    description=(
        "dem4water estimate a law between altitude and "
        "surface to estimate water volume in dam."
    ),
    entry_points={
        "console_scripts": [
            "dem4water=dem4water.cli:main",
            "config_template=dem4water.tools.generate_default_configuration:main",

        ],
    },
    install_requires=requirements,
    long_description=readme,
    include_package_data=True,
    keywords="dem4water",
    name="dem4water",
    packages=find_packages(include=["dem4water", "dem4water.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://gitlab.cnes.fr/si2a/dem4water",
    version="0.1.0",
    zip_safe=False,
)
