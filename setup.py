import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hivdrm",
    version="0.0.1",
    author="Sergey Naumenko",
    description="Detect HIV Drug Resitant Mutations using amplicon sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    project_urls={
        "Bug Tracker": "https://github.com/bcbio/hivdrm/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "hivdrm"},
    packages=setuptools.find_packages(where="hivdrm"),
    python_requires=">=3.6",
)