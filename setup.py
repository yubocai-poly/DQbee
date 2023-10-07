from setuptools import setup, find_packages
import os

VERSION = '0.1.5'
DESCRIPTION = 'DQode is a Python Library to compute the optimal inner-quadratic quadratization and dissipative quadratization of a given polynomial ODE system.'

setup(
    name="DQbee",
    version=VERSION,
    author="Yubo Cai, Gleb Pogudin",
    author_email="yubocai0811@gmail.com, pogudin.gleb@gmail.com",
    maintainer="Yubo Cai",
    maintainer_email="yubocai0811@gmail.com",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=open('README.md','r',encoding="UTF8").read(),
    packages=find_packages(),
    install_requires=['sympy', 
                      'numpy', 
                      'scipy', 
                      'matplotlib', 
                      'tbcontrol'],
    keywords=['quadratization', 'differential equation', 'symbolic computing'],
    url="https://github.com/yubocai-poly/Dissipative-Quadratiation-Package",
    project_urls={
        "Code": "https://github.com/yubocai-poly/Dissipative-Quadratiation-Package/",
        "Issue tracker": "https://github.com/yubocai-poly/Dissipative-Quadratiation-Package/issues",
    },
    license="MIT",
    classifiers= [
        "Operating System :: OS Independent",
        "Topic :: Text Processing :: Indexing",
        "Topic :: Utilities",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
)