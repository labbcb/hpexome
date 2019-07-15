import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='hpexome',
    version='1.0.0',
    author="Welliton Souza",
    author_email="well309@gmail.com",
    description="An automated workflow for processing whole-exome sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/labbcb/hpexome",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'Click'
    ],
    entry_points='''
        [console_scripts]
        hpexome=hpexome.hpexome:hpexome
    '''
)
