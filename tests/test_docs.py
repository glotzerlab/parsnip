# # import subprocess

# # import pytest


# # @pytest.mark.sphinx()
# # def test_sphinx_doctests():
# #     """Run Sphinx's doctests."""
# #     result = subprocess.run(
# #         ["sphinx-build", "-b", "doctest", "doc/source/", "_build/doctest"],
# #         capture_output=True,
# #         text=True,
# #     )
# #     print(result.stdout)
# #     print(result.stderr)

# #     assert result.returncode == 0, "Sphinx doctests failed"
# import pytest
# from sphinx.application import Sphinx

# @pytest.mark.sphinx
# def test_sphinx_doctests():
#     """Run Sphinx doctests programmatically."""
#     src_dir = "doc/source/"  # Path to the directory containing your .rst files
#     build_dir = "_build"  # Path to store the build output
#     doctest_dir = f"{build_dir}/doctest"

#     # Create a Sphinx application instance
#     app = Sphinx(
#         srcdir=src_dir,
#         confdir=src_dir,
#         outdir=doctest_dir,
#         doctreedir=f"{build_dir}/doctrees",
#         buildername="doctest",
#     )

#     # Build the documentation and run doctests
#     app.build()

#     # Check for build errors
#     if app.statuscode != 0:
#         pytest.fail(f"Sphinx doctests failed with status code {app.statuscode}")
