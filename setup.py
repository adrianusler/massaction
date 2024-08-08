from massaction import __version__
import setuptools

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="massaction",
    version=__version__,
    description="A python package for the numerical solution of mass-action laws with constraints for the description of chemical reactions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adrianusler/massaction",
    author="Adrian L. Usler",
    author_email="adrian.usler@rwth-aachen.de",
    license="MIT License",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "scipy",
    ],
    python_requires="~=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
