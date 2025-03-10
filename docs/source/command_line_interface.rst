Command Line Interface
#################

The **Command Line Interface (CLI)** can be used to run analysis from the command line or a Python console.

.. contents:: Contents of this Page
   :local:


Run a Simulation
****************

To run a simulation, we first need to create a project.

Create New Project
^^^^^^^^^^^^^^^^^^

First, the relevant modules should be imported.
To create a new project,

.. code-block:: python

   import numpy as np
   from cli.cavity import Cavity, Cavities
   import matplotlib.pyplot as plt
   plt.ion()


The geometric inputs are lists or numpy array. Therefore ``numpy`` module should be imported.
``matplotlib`` is and interactive mode enabled for plotting.

.. warning::
    Interactive plotting must be enabled using ``plt.ion()`` otherwise the plots are not displayed.

The cavity can now be defined as

.. code-block:: python

   mid_cell = np.array([42, 42, 12, 19, 35, 57.7, 103.3])
   cav = Cavity(2, mid_cell)
   cav.save(files_path='D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\ConsoleTest')

.. warning::
   It is important to save the cavity first before proceeding to run analysis on the cavity

The example above defines a 2-cell cavity

Eigenmode Analysis
^^^^^^^^^^^^^^^^^^

To run the eigenmode analysis use

.. code-block:: python

   cav.run_eigenmode()

This uses default settings for the analysis.


Tune
^^^^

To tune the cavity, use

.. code-block:: python

   <cavity_name>.run_tune('Req', freq=1300)

This uses default settings for the analysis.

Wakefield Analysis
^^^^^^^^^^^^^^^^^^^^

To run wakefield analysis

.. code-block:: python

   <cavity_name>.run_wakefield()

To run the wakefield impedance analysis for a particular working point, the working points has to be defined first as
a dictionary as shown below

.. code-block:: json

   wp = {
            'H': {
                'I0 [mA]': '26.7',
                'Nb [1e11]': '1.51',
                'sigma_z (SR/BS) [mm]': '2.5/4.45'
            }
   }

Plotting
********

.. code-block:: python

   cavs.plot_compare_fm_bar()
