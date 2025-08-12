.. _install:

=========================
Installation Instructions
=========================

DRAIN-LID is a MATLAB-based model (tested with MATLAB R2021a and newer).
No Python installation is required to run the model itself.
Python is only needed if you want to rebuild this documentation or use optional post-processing scripts.

MATLAB Installation
-------------------
1. **Clone (or download) the repository**:
   ::
      git clone https://github.com/marcusnobrega-eng/LID_Tool.git

2. **Add DRAIN-LID to your MATLAB path**:
   - In MATLAB:
   ::
        addpath(genpath('path/to/LID_Tool'))

3. **Test your installation** by running the quickstart example:
   - In MATLAB:
   ::
        addpath(genpath('path/to/LID_Tool/Model_Configurations')
   ::
        run('test_model.m')

If plots are generated and saved in `./Outputs/Example01`, your installation is working.


Developer Notes
---------------
- Test cases are located in the `examples/` folder.
- See the :doc:`../quickstart/index` page for a minimal run example.
