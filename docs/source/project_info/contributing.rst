Contributing
=======================================

We welcome and appreciate any and all contributions.
Any effort will be given due credit.


Feedback
----------------------
To report a bug or propose a new feature, please open an issue `here <https://github.com/PNNL-m-q/mzapy/issues>`_. 

When reporting bugs, please include:

* Your operating system and version.
* The version(s) of the package dependency (or dependencies) that produce the bug
* Detailed steps to reproduce the bug
* Any details about your local setup that might be helpful in troubleshooting

When proposing a new feature:

* Explain in detail how it would work
* Keep the scope as narrow as possible
* Where possible, include code/pseudocode snippets to explain implementation or expected behavior


Fixing Bugs or Implementing Features
------------------------------------------

Look through the Github issues for bugs or proposed features.
If you want to start working on a bug or feature, please write short message on the issue tracker 
to prevent duplicate work. 
Bugfixes and implemented features should be submitted as pull requests.


Writing Documentation
-------------------------------
We welcome the addition of more documentation, or improvements to existing documentation. 
This include the documentation you are reading now, docstrings, tutorials, or even on the web as blog posts or articles.
Please reference the documentation guidelines below.
Changes to documentation should be submitted as pull requests.


Submitting a Pull Request
-------------------------

To fix bugs, add new features, or update documentation you need to create a pull request.
Any changes you make to your local copy of the code will be made easily available for review and 
eventual integration into the code base.
To create a pull request:

#. Create a GitHub account
#. Fork the repository
#. Clone your fork locally
#. Go to the created ``mzapy`` folder with :code:`cd mzapy`
#. Create a new branch (based on current ``dev`` branch) with: :code:`git checkout -b <your_github_username>_<type>_<descriptive_branch_name> dev`, where ``<type>`` is "bugfix", "feature", or "docs"
#. Make your changes to the code or documentation
#. Run :code:`git add --all` to add all the changed files to the commit
#. To commit the added files, run :code:`git commit -m '<text>'` or :code:`git commit`. The latter will open a command line editor to write a commit message. Commit messages should have a descriptive line header (<50 characters, including spaces), followed by an empty line, and then a description of what you did and why.
#. Push your changes to your ``mzapy`` fork by running :code:`git push origin <your_github_username>_<type>_<descriptive_branch_name>`
#. If you then navigate to the webpage for your mzapy fork, you should see a "Pull request" link in the sidebar. Choose the relevant pull request from the menu and click the "Create Pull Request" button. Ensure the pull request target branch is `<your_github_username>_<type>_<descriptive_branch_name>!`

If you want to create more pull requests, first run :code:`git checkout dev` and then start at step 5 with a new branch name.


Versioning Scheme
----------------------------
The version of ``mzapy`` is stored in the ``__version__`` variable in ``mzapy/__init__.py``.
The versioning scheme follows the format *mza_version.major_version.minor_version*:

* *mza version* - this is kept in lock-step with the version of the underlying MZA format
* *major version* - increments after every time changes in the ``dev`` branch are merged into the ``main`` branch
* *minor version* - increments after every time changes in individual bugfix/feature/docs branches are merged into the ``dev`` branch

.. note::

    When creating a branch to implement a bugfix, new feature, or add documentaion, append ".<your_github_username>0" 
    to the current ``__version__`` variable in ``mzapy/__init__.py``, then increment the number at the end with each 
    commit you make to that branch.
    

Coding Style
-----------------------------

Naming Conventions
******************************
Functions and classes should be named in a way that describes what they do and whether they are internal 
(not intended to be part of the public-facing API) or external (intended to be part of the public-facing API). 
Internal function/class names should be prepended with "_". Example:

.. code-block:: python3

    # this is an internal function, not meant to be part of the public-facing API
    def _add(a, b):
        """ adds two integers """
        return a + b

    # this is an external function, meant to be part of the public-facing API
    def sum_pairwise(x, y):
        """ 
        returns the pairwise sums of integers from two lists

        Parameters
        ----------
        x : list(int)
        y : list(int)
            lists of integers, must be same length

        Returns
        -------
        sums : list(int)
            list of pairwise sums
        """
        # uses the internal _add function
        return [_add(a, b) for a, b in zip(x, y)]


Docstring Format
****************************
Detailed docstrings must be included in all functions/classes (both internal and external) in ``mzapy``. Docstring 
format loosely follows the `numpydoc style <https://numpydoc.readthedocs.io/en/latest/format.html>`_, refer to 
existing docstrings for specific examples. Generally, all functions should include a description and parameters/returns 
sections (if applicable) as in the following example:

.. code-block:: python3
    
    def foo(a, b, c=None, d=1234):
        """
        Give a brief description about what the function does, what inputs it takes, and what outputs it produces
        
        Parameters
        ----------
        a : int
            parameter a description 
        b : float
            parameter b description
        c : str, optional
            parameter c description, indicate behaviors when c parameter is provided/not provided
        d : int, default=1234
            parameter d description, if the default value has some significance describe that here

        Returns
        -------
        x : int
            description of return value, add more entries if the function returns more than one thing
        """
        ...


.. note::

    The Parameters and Returns sections may be omitted if a function does not take parameters and/or produce a 
    return value. They may also be omitted if the function performs a trivial enough task that parameters and return
    values can easily be inferred from the description.


Adding Entries to Sphinx Documentation
**********************************************
All public-facing functions/classes should have entries in the Sphinx documentation source files in order for their 
docstrings to be incorporated into the HTML documentation. For example, if the ``sum_pairwise`` function in the example 
above were implemented in the ``mzalib/isotopes.py`` module, then the following entry should be added to the 
appropriate section in the ``docs/source/isotopes.rst`` documentation source file:

.. code-block::

    Module Reference
    ------------------------------
    
    .. autofunction :: mzapy.isotopes.sum_pairwise



