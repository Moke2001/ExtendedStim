from pathlib import Path
from setuptools import setup, find_packages

HERE = Path(__file__).parent
README = (HERE / "README.md").read_text(encoding="utf-8")

setup(
    name="extendedstim",
    version="0.1.0",
    description="ExtendedStim: tools for testing quantum circuits, quantum error-correction codes.",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Moke",
    author_email="moke2001@whu.edu.cn",
    license="MIT",
    packages=find_packages(exclude=("tests", "docs")),
    include_package_data=True,
    python_requires=">=3.12",
    install_requires=["qiskit", "numpy", "stim", "scipy", "matplotlib","stimbposd","pylatexenc","qutip","galois"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

