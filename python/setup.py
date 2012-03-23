import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "biweb",
    version = "1.0.0",
    author = "Timothee Poisot",
    author_email = "timothee.poisot@uqar.ca",
    description = ("A Python module to work on bipartite networks"
                   "of ecological interactions."),
    license = "GNU GPL",
    keywords = "ecology bipartite networks bioinformatics",
    url = "http://tpoisot.github.com/biweb/",
    packages=['biweb','biweb.base','biweb.dataIO','biweb.null','biweb.nes','biweb.mod','biweb.graphs','biweb.tests','biweb.spe'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Intended Audience :: Science/Research",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization"
        ],
    )
