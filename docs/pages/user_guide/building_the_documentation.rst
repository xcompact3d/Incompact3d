Building the documentation
==========================

- After cloning the Xcompact3d's repository to your system, navigate to its folder;
- Make sure to install all requirements with:

  .. code-block:: bash

     pip install -r ./docs/requirements.txt

- It is time to build it:

  .. code-block:: bash

     sphinx-autobuild docs docs/_build/html

  Following `sphinx-autobuild's instructions <https://pypi.org/project/sphinx-autobuild/>`_, This will start a server at http://127.0.0.1:8000 and start watching for changes in the ``docs/`` directory. When a change is detected in ``docs/``, the documentation is rebuilt and any open browser windows are reloaded automatically. ``KeyboardInterrupt`` (ctrl+c) will stop the server.
