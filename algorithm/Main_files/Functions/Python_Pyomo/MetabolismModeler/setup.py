from setuptools import setup, find_packages

setup(
    name="myPackage",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        # add your dependencies here
    ],
    include_package_data=True,
    description="A description of your package",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://your.package.url",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
