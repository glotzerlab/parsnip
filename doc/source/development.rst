.. _development:

=================
Development Guide
=================


All contributions to **parsnip** are welcome!
Developers are invited to contribute to the framework by pull request to the package repository on `GitHub`_, and all users are welcome to provide contributions in the form of **user feedback** and **bug reports**.
We recommend discussing new features in form of a proposal on the issue tracker for the appropriate project prior to development.

.. _github: https://github.com/glotzerlab/parsnip

General Guidelines
==================

All code contributed to **parsnip** must adhere to the following guidelines:

  * Use a two branch model of development:

    - Most new features and bug fixes should be developed in branches based on ``main``.
    - API incompatible changes and those that significantly change existing functionality should be based on ``breaking``
  * Hard dependencies (those that end users must install to use **parsnip**) are *strongly* discouraged, and should be avoided where possible. Numpy is the sole exception to this, as it is already included in most scientific computing software stacks.
  * Additional dependencies required by developers (those used to run tests or build docs) are undesirable but allowed.
  * All code should adhere to the source code conventions and satisfy the documentation and testing requirements discussed below.

As portability is a primary feature of **parsnip**, tests are run run on Python versions 3.9 and later. However, first class support should only be expected for versions covered by `NEP 29`_.

.. _NEP 29: https://numpy.org/neps/nep-0029-deprecation_policy.html


Style Guidelines
----------------

The **parsnip** package adheres to a reasonably strict set of style guidelines.
All code in **parsnip** should be formatted using `ruff`_ via pre-commit. This provides an easy workflow to enforce a number of style and syntax rules that have been configured for the project.

.. tip::

    `pre-commit`_ has been configured to run a number of linting and formatting tools. It is recommended to set up pre-commit to run automatically:

    .. code-block:: bash

        python -m pip install pre-commit
        pre-commit install # Set up git hook scripts

    Alternatively, the tools can be run manually with the following command:

    .. code-block:: bash

        git add .; pre-commit run

.. _ruff: https://docs.astral.sh/ruff/
.. _pre-commit: https://pre-commit.com/


Documentation
-------------

API documentation should be written as part of the docstrings of the package in the `Google style <https://google.github.io/styleguide/pyguide.html#383-functions-and-methods>`__.

Docstrings are automatically validated using `pydocstyle <http://www.pydocstyle.org/>`_ whenever the ruff pre-commit hooks are run.
The `official documentation <https://parsnip.readthedocs.io/>`_ is generated from the docstrings using `Sphinx <http://www.sphinx-doc.org/en/stable/index.html>`_.

In addition to API documentation, inline comments are strongly encouraged.
Code should be written as transparently as possible, so the primary goal of documentation should be explaining the algorithms or mathematical concepts underlying the code.
Avoid comments that simply restate the nature of lines of code: for example, the comment "set up regex pattern" is uninformative, since the code itself should make this obvious, *e.g*, ``re.compile(r"^(_\w+)\s+(\d+)")``.
On the other hand, the comment "read an underscore-prefixed word and a numeric value seperated by whitespace" is instructive.
Multiline comments for regex strings may sometimes be necessary.


Unit Tests
----------

All code should include a set of tests which test for correct behavior.
All tests should be placed in the ``tests`` folder at the root of the project.
In general, most parts of parsnip primarily require `unit tests <https://en.wikipedia.org/wiki/Unit_testing>`_, but where appropriate `integration tests <https://en.wikipedia.org/wiki/Integration_testing>`_ are also welcome. Core functions should be tested against the sample CIF files included in ``tests/sample_data``.
Tests in **parsnip** use the `pytest <https://docs.pytest.org/>`__ testing framework.
To run the tests, simply execute ``pytest`` at the root of the repository.
