from setuptools import setup

setup(
    name="PanR2",
    version="0.1.2",
    author="Tasnimul Arabi Anik",
    author_email="arabianik987@gmail.com",
    description="A Python tool for panresistome analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Tasnimul-Arabi-Anik/PanR2",
    scripts=["bin/panr"],  # Include the script from the bin directory
    install_requires=[
	"pandas>=1.3",
	"matplotlib>=3.5",
	"seaborn>=0.11",
	"numpy>=1.21",
	"scipy>=1.7",
	"plotly>=5.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
