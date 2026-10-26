.. _vs_app_apk:

=============================
Visual Sudoku App for Android
=============================

.. include:: menu_backto.rst

A visual sudoku solver based on cost function networks
------------------------------------------------------

  This application solves the sudoku problem from a smartphone by reading the
  grid using its camera. The cost function network solver toulbar2 is used to
  deal with the uncertainty on the digit recognition produced by the neural
  network. This uncertainty, combined with the sudoku logical rules, makes it
  possible to correct perceptual errors. It is particularly useful in the case
  of hand-written digits or poor image quality. It is also possible to solve a
  partially filled-in grid with printed and hand-written digits. The solver
  will always suggest a valid solution that best adapts to the retrieved digit
  information. It will naturally detect (a small number of) errors in a
  partially filled-in grid and could be used later as a diagnosis tool (future
  work). This software demonstration emphasizes the tight relation between
  constraint programming, computer vision, and deep learning.

  We used the open-source C++ solver
  `toulbar2 <https://github.com/toulbar2/toulbar2>`_ in order to find the
  maximum a posteriori solution of a constrained probabilistic graphical model.
  With its dedicated numerical (soft) local consistency bounds, toulbar2
  outperforms traditional CP solvers on this problem. Grid perception and cell
  extraction are performed by the computer vision library
  `OpenCV <https://opencv.org>`_. Digit recognition is done by **Keras** and
  **TensorFlow**. The current android application is written in Python using
  the `Kivy <https://kivy.org>`_ framework. It is inspired from a
  `tutorial <https://pyimagesearch.com/2020/08/10/opencv-sudoku-solver-and-ocr>`_
  by Adrian Rosebrock. It uses the `ws <http://147.100.179.250>`_ RESTful web
  services in order to run the solver.  

  See also : :ref:`vs_app`.

Source Code
-----------

  `GitHub code <https://github.com/toulbar2/visualsudoku>`_
  |github_logo_visualsudokuapk|

Download and Install
--------------------

  To install the 'Visual Sudoku' application on smartphone :

   1) **Download** the **visualsudoku-release.apk** APK file from Github
      repository : 

       |qrcode_visualsudokuapk| 
       `https://github.com/toulbar2/visualsudoku/releases/latest <https://github.com/toulbar2/visualsudoku/releases/latest>`_

   2) Click on the downloaded **visualsudoku-release.apk** APK file to ask
      for **installation** 
      (*you have to accept to 'install anyway' from unknown developer*).

   3) In your parameter settings for the app, give permissions to the
      'Visual Sudoku' application
      (smartphone menu 'Parameters' > 'Applications' > 'Visual Sudoku') :
      allow camera (required to capture grids), files and multimedia contents
      (required to save images as files). Re-run the app.

  Warnings :

    - The application may fail at first start and you may have to launch it
      twice.
    - While setting up successfully, the application should have created
      itself the required 'VisualSudoku' folder (under the smartphone
      'Internal storage' folder) but if not, you will have to create it
      by yourself manually.
    - Since the application calls a web service, an internet connection is required.

Description
-----------

  The 'SETTINGS' menu allows to save grids or solutions as image files ('savinginputfile', 'savingoutputfile' parameters) and to access to some 'expert' parameters in order to enhance the resolution process ('keep', 'border', 'time' parameters). 

  .. include:: vsapp/vsappapk_overview.rst


.. |github_logo_visualsudokuapk| image:: /_static/img/logo-github.png
   :width: 30
   :alt: https://github.com/toulbar2/visualsudoku
   :target: https://github.com/toulbar2/visualsudoku

.. |qrcode_visualsudokuapk| image:: /_static/img/qr-code_visualsudokuapk.png
   :width: 100
   :alt: `url_download_visualsudokuapk`_
   :target: `url_download_visualsudokuapk`_

.. _url_download_visualsudokuapk: https://github.com/toulbar2/visualsudoku/releases/latest

