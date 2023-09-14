from setuptools import find_packages, setup

setup(
    name='valueengineering',
    packages=find_packages(include=["valueengineering"]),
    version='0.1.0',
    description='Tools developed by Rasmus Rune Pedersen at Value Engineering Aps',
    author='Me',
    license='MIT',
    install_requires=['numpy']
)
