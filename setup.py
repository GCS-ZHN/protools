#!/usr/bin/env python
from setuptools import find_packages, setup
import versioneer

PACKAGE_DIR = "."

setup(
    package_dir = {"": PACKAGE_DIR},
    packages = find_packages(where=PACKAGE_DIR, exclude=['test', 'test.*']),
    include_package_data=True,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)
