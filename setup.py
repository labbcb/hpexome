import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='HPexome',
    version='1.1.0',
    author="Welliton Souza",
    author_email="well309@gmail.com",
    description="An automated tool for processing whole-exome sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://bcblab.org/hpexome",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3"
    ],
    install_requires=[
        'Click'
    ],
    entry_points='''
        [console_scripts]
        hpexome=hpexome.hpexome:hpexome
    ''',
    project_urls={
        "Source Code": "https://github.com/labbcb/hpexome",
        "Bug Tracker": "https://github.com/labbcb/hpexome/issues"
    }
)
