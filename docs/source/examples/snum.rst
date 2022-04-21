.. _snum:

==================
Sudoku Application
==================

.. include:: menu_backto.rst

Brief description
=================

A Sudoku code returning a sudoku partial grid (sudoku problem) and the corresponding completed grid (sudoku solution), such as :download:`partial and completed grids<snum/mygrids.txt>`.

The verbose version, that further gives a detailed description of what the program does, could be useful as tutorial example. Example : 
:download:`partial and completed grids with explanations<snum/mygrids_details.txt>`.

Available
=========

- **As a web service** :

  You can run the software directly from your browser as a web service :

  Grids information is returned into the output stream.
  The **returned_type** parameter of the web service allows to choose
  how to receive it :

    - *returned_type=stdout.txt* : to get the output stream as a .txt file.
    - *returned_type=run.zip* : to get the .zip run folder containing the output stream __WS__stdout.txt (*+ the error stream __WS__stderr.txt that may be useful to investigate*).

  - web service to get one sudoku grids (both partial and completed) :
    `api/ui/sudoku <http://147.100.179.250/api/ui/sudoku>`_
    |qrcode_ui_sudoku|

  - web service to further get a detailed description
    of what the program does :
    `api/ui/sudoku/tut (verbose version) <http://147.100.179.250/api/ui/sudoku/tut>`_
    |qrcode_ui_sudokutut| 

  The **sudoku web services**, hosted by the `ws <http://147.100.179.250>`_ web services (based on HTTP protocol), can be called by many other ways : from a browser (like above), from any softwares written in a language supporting HTTP protocol (Python, R, C++, Java, Php...), from command line tools (cURL...)...

  For example, calling the sudoku web services from a terminal by cURL :

  - Commands (*replace indice value by any value in 1...17999*) : ::

      curl --output mygrids.txt -F 'indice=778' -F 'returned_type=stdout.txt' http://147.100.179.250/api/tool/sudoku

      curl --output myrun.zip -F 'indice=778' -F 'returned_type=run.zip' http://147.100.179.250/api/tool/sudoku

      # verbose version

      curl --output mygrids_details.txt -F 'indice=778' -F 'returned_type=stdout.txt' http://147.100.179.250/api/tool/sudoku/tut

      curl --output myrun_details.zip -F 'indice=778' -F 'returned_type=run.zip' http://147.100.179.250/api/tool/sudoku/tut

  - Responses corresponding with the requests above :

    - :download:`mygrids.txt<snum/mygrids.txt>`
    - __WS__stdout.txt into myrun.zip has the same content as
      :download:`mygrids.txt<snum/mygrids.txt>`
    - :download:`mygrids_details.txt<snum/mygrids_details.txt>` (*__WS__stdout.txt into myrun_details.zip has the same content*)
    - __WS__stdout.txt into myrun_details.zip has the same content as
      :download:`mygrids_details.txt<snum/mygrids_details.txt>`

.. |qrcode_ui_sudoku| image:: /_static/img/qr-code_ui-sudoku.png
   :width: 100
   :alt: `url_ui_sudoku`_
   :target: `url_ui_sudoku`_

.. _url_ui_sudoku: http://147.100.179.250/api/ui/sudoku

.. |qrcode_ui_sudokutut| image:: /_static/img/qr-code_ui-sudoku-tut.png
   :width: 100
   :alt: `url_ui_sudokutut`_
   :target: `url_ui_sudokutut`_

.. _url_ui_sudokutut: http://147.100.179.250/api/ui/sudoku/tut

