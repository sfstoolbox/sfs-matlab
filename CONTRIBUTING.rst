Contributing
------------

If you find errors, omissions, inconsistencies or other things that need
improvement, please create an issue or a pull request at
https://github.com/sfstoolbox/sfs-matlab/.
Contributions are always welcome!

Building the Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you make changes to the documentation, you can re-create the HTML pages
using Sphinx_.
You can install it and a few other necessary packages with::

   python3 -m pip install -r doc/requirements.txt --user

To create the HTML pages, use::

   sphinx-build -b html -d ./_build/doctrees . ./_build/html-preview/

The generated files will be available in the directory
``build/sphinx/html-preview/``.

.. _Sphinx: http://sphinx-doc.org/

Creating a New Release
^^^^^^^^^^^^^^^^^^^^^^

New releases are made using the following steps:

#. Bump version number in ``SFS_version.m``
#. Update ``NEWS``
#. Commit those changes as "Release x.y.z"
#. Create an (annotated) tag with ``git tag -a x.y.z``
#. Check that all author details in ``.zenodo.json`` are correct
#. Push the commit and the tag to Github and `add release notes`_ containing a
   link to the documentation with https://sfs-matlab.readthedocs.io/en/x.y.z and
   the bullet points from ``NEWS``
#. Check that the new release was built correctly on RTD_, delete the "stable"
   version and select the new release as default version
#. Check that the new release was registered at Zenodo_ and edit the release
   notes on Github to include the DOI badge

.. _add release notes: https://github.com/sfstoolbox/sfs-matlab/tags
.. _RTD: https://readthedocs.org/projects/sfs-matlab/builds/
.. _Zenodo: https://zenodo.org/search?page=1&size=20&q=SFS%20Toolbox&sort=bestmatch
