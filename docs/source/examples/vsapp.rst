.. _vs_app:

=========================
Visual Sudoku Application
=========================

.. include:: menu_backto.rst

Brief description
=================

An automatic Sudoku puzzle **solver** using **OpenCV**, **Deep Learning**, and **Optical Character Recognition** (**OCR**).

Available
=========

- **Software** adapted by Simon de Givry (@ INRAE, 2022)
  in order to use **toulbar2** solver, from a
  `tutorial <https://pyimagesearch.com/2020/08/10/opencv-sudoku-solver-and-ocr>`_
  by Adrian Rosebrock (@ PyImageSearch, 2022) :

  `GitHub code <https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/tree/master/crafted/visualsudoku>`_
  |github_logo_vsudoku|
  
- **As a web service** :

  - You can run this software directly from your browser as a
    `web service <http://147.100.179.250/api/ui/vsudoku>`_
    |qrcode_ui_vsudoku| 
   
  - or from your smartphone (install the following android application ; in your parameter settings for the app, allows accessing to camera and files ; rerun the app) `Visual Sudoku App <http://147.100.179.250/visualsudoku-release.apk>`_

  - Other ways to call the web service :

    The **visual sudoku web service**, hosted by the `ws <http://147.100.179.250>`_ web services (based on HTTP protocol), can be called by many ways : from a browser (like above), from any softwares written in a language supporting HTTP protocol (Python, R, C++, Java, Php...), from command line tools (cURL...)...

    Example from a terminal by cURL (*replace mygridfilename.jpg by your own image file name*) : ::

      curl --output mysolutionfilename.jpg -F 'file=@mygridfilename.jpg' -F 'keep=40' -F 'border=15' http://147.100.179.250/api/tool/vsudoku

.. |github_logo_vsudoku| image:: /_static/img/logo-github.png
   :width: 30
   :alt: https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/tree/master/crafted/visualsudoku
   :target: https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/tree/master/crafted/visualsudoku

.. |qrcode_ui_vsudoku| image:: /_static/img/qr-code_ui-vsudoku.png
   :width: 100
   :alt: `url_ui_vsudoku`_
   :target: `url_ui_vsudoku`_

.. _url_ui_vsudoku: http://147.100.179.250/api/ui/vsudoku

