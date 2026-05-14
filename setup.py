from setuptools import setup, find_packages

setup(
    name="Stellar-Interior-Model",  # El nombre que tendrá tu paquete
    version="1.0.0",
    author="Pablo Paez Ramos",
    author_email="pablopaezramos2003@gmail.com",
    description="Modeling the stellar interior of a star with a convective core and a radiative envelope based on its chemical composition and mass.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/pablo-paez/Stellar-Interior-Model",
    packages=find_packages(), # Busca automáticamente todas las carpetas con __init__.py
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Astrophysics", # Cambia según tu TFG
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "requests",
        "matplotlib"
    ],
)