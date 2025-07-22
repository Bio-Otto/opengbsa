Contributing to OpenGBSA
=======================

We welcome contributions from the community! You can help by reporting issues, requesting features, improving documentation, or contributing code.

How to Report Issues or Request Features
----------------------------------------
.. tip::
   Open an issue on GitHub: https://github.com/Bio-Otto/opengbsa/issues
   
   Clearly describe the problem or feature request, and include error messages, logs, and environment details if relevant.

How to Contribute Code
----------------------
1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   .. code-block:: bash

      git clone https://github.com/your-username/opengbsa.git
      cd opengbsa

3. **Create a new branch** for your feature or fix:
   .. code-block:: bash

      git checkout -b my-feature

4. **Make your changes**
   - Follow PEP8 coding style
   - Use clear, descriptive commit messages
   - Add or update docstrings for all public functions/classes

5. **Run tests** to ensure nothing is broken:
   .. code-block:: bash

      pytest test/

6. **Build the documentation** (optional, for doc changes):
   .. code-block:: bash

      cd docs
      make html

7. **Push your branch** to your fork:
   .. code-block:: bash

      git push origin my-feature

8. **Open a pull request** on GitHub
   - Describe your changes and reference any related issues

.. note::
   Keep pull requests focused and small if possible. This makes review easier and faster.

How to Contribute Documentation
-------------------------------
.. tip::
   Edit or add files in the ``docs/`` directory (reStructuredText or Markdown). Build the docs locally to check formatting, and follow the style of existing documentation.

Code Style and Best Practices
-----------------------------
.. important::
   - Follow PEP8 (use tools like flake8 or black)
   - Write clear docstrings (Google or NumPy style)
   - Add tests for new features or bugfixes

Thank you for helping improve OpenGBSA! 