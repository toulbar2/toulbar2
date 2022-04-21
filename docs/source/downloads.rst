.. _downloads:

=========
Downloads
=========

Packages
--------

- to install **toulbar2** using the package manager in Debian and Debian derived Linux distributions (Ubuntu, Mint,...): ::

    apt install toulbar2

Binaries
--------

- **Latest release toulbar2 binaries**

  `Linux 64bit <https://github.com/toulbar2/toulbar2/releases/download/v1.1.1/toulbar2>`_ |
  `MacOs 64bit <https://github.com/toulbar2/toulbar2/releases/download/v1.1.1/toulbar2mac>`_ |
  `Windows 64bit <https://github.com/toulbar2/toulbar2/releases/download/v1.1.1/toulbar2.exe>`_

Python package
--------------

- **pytoulbar2** module for Linux and MacOS on `PyPI <https://pypi.org/project/pytoulbar2>`_ 

Open-source code
----------------

- `toulbar2 on GitHub <https://github.com/toulbar2/toulbar2>`_

Docker images
-------------

- In `Toulbar2 Packages <https://github.com/toulbar2?tab=packages&repo_name=toulbar2>`_ :

  - **toulbar2** : 
    Docker image containing toulbar2 and its pytoulbar2 Python API
    *(installed from sources with cmake options -DPYTB2=ON and -DXML=ON)*.
    Install from the command line: ::

      docker pull ghcr.io/toulbar2/toulbar2/toulbar2:master

  - **pytoulbar2** :
    Docker image containing pytoulbar2 the Python API of toulbar2
    *(installed with python3 -m pip install pytoulbar2)*.
    Install from the command line: ::

      docker pull ghcr.io/toulbar2/toulbar2/pytoulbar2:master

