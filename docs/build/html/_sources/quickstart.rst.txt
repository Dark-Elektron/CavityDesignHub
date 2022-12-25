Quick Start Guide
#################

This **Quick Start Guide** gives a fast and simple introduction into PoC. All
topics can be found in the :ref:`Using PoC <USING>` section with much more
details and examples.

.. contents:: Contents of this Page
   :local:


.. _QUICK:Requirements:

Requirements and Dependencies
*****************************

The CaDH is a GUI for coordinating different electromagnetic codes.
It also comes with some scripts to ease most of the common tasks. CaDH uses
Python 3 as a platform independent scripting environment.
All Python scripts are wrapped in Bash or PowerShell scripts,
to hide some platform specifics of Darwin, Linux or Windows.
See :ref:`USING:Require` for further details.


.. rubric:: CaDH requires:

* The **Python 3** programming language and runtime, if you want to use PoC's infrastructure.
* A shell to execute shell scripts:

  * **Bash** on Linux and OS X
  * **PowerShell** on Windows


.. rubric:: CaDH optionally requires:

* **Git** command line tools or
* **Git User Interface**, if you want to check out the latest 'master' or 'release' branch.

All dependencies are contained in requirements.txt


.. _QUICK:Download:

Download
********

The CaDH-Library can be downloaded as a `zip-file <https://github.com/VLSI-EDA/PoC/archive/master.zip>`_
(latest 'master' branch), cloned with ``git clone`` or embedded with
``git submodule add`` from GitHub. GitHub offers HTTPS and SSH as transfer
protocols. See the :ref:`Download <USING:Download>` page for further
details. The installation directory is referred to as ``CaDHRoot``.

+----------+---------------------------------------------------------------------+
| Protocol | Git Clone Command                                                   |
+==========+=====================================================================+
| HTTPS    | |
+----------+---------------------------------------------------------------------+
| SSH      |  |
+----------+---------------------------------------------------------------------+


.. _QUICK:Configuration:

Configuring CaDH  on a Local System
*********************************

To explore CaDH's full potential, it's required to configure some paths and
synthesis or simulation tool chains. The following commands start a guided
configuration process. Please follow the instructions on screen. It's possible
to relaunch the process at any time, for example to register new tools or to
update tool versions. See :ref:`Configuration <USING:PoCConfig>` for
more details. Run the following command line instructions to configure PoC on
your local system:

.. code-block:: PowerShell

   cd CaDHRoot
   .\poc.ps1 configure


.. _QUICK:Integration:

Integration
***********

The CaDH-Library is meant to be integrated into other HDL projects. Therefore
it's recommended to create a library folder and add the CaDH-Library as a Git
submodule. After the repository linking is done, some short configuration
steps are required to setup paths, tool chains and the target platform. The
following command line instructions show a short example on how to integrate
PoC.

.. rubric:: 1. Adding the Library as a Git submodule

The following command line instructions will create the folder ``lib\PoC\`` and
clone the PoC-Library as a Git `submodule <http://git-scm.com/book/en/v2/Git-Tools-Submodules>`_
into that folder. ``ProjectRoot`` is the directory of the hosting Git. A detailed
list of steps can be found at :doc:`Integration </UsingPoC/Integration>`.

.. code-block:: PowerShell

   cd ProjectRoot
   mkdir lib | cd
   git submodule add https://github.com/VLSI-EDA/PoC.git PoC
   cd PoC
   git remote rename origin github
   cd ..\..
   git add .gitmodules lib\PoC
   git commit -m "Added new git submodule PoC in 'lib\PoC' (PoC-Library)."


.. _QUICK:RunSimulation:

Run a Simulation
****************

The following quick example uses the GHDL Simulator to analyze, elaborate and
simulate a testbench for the module ``arith_prng`` (Pseudo Random Number
Generator - PRNG). The VHDL file ``arith_prng.vhdl`` is located at
``PoCRoot\src\arith`` and virtually a member in the `PoC.arith` namespace.
So the module can be identified by an unique name: ``PoC.arith.prng``, which is
passed to the frontend script.

.. rubric:: Example:

.. code-block:: PowerShell

   cd CaDHRoot
   python main_control.py

The CLI command ``ghdl`` chooses *GHDL Simulator* as the simulator and
passes the fully qualified PoC entity name ``PoC.arith.prng`` as a parameter
to the tool. All required source file are gathered and compiled to an
executable. Afterwards this executable is launched in CLI mode and its outputs
are displayed in console:

.. image:: /_static/images/ghdl/arith_prng_tb.posh.png
   :target: /_static/images/ghdl/arith_prng_tb.posh.png
   :alt: PowerShell console output after running PoC.arith.prng with GHDL.

Each testbench uses PoC's simulation helper packages to count asserts and to
track active stimuli and checker processes. After a completed simulation run,
an report is written to STDOUT or the simulator's console. Note the line
``SIMULATION RESULT = PASSED``. For each simulated PoC entity, a line in the
overall report is created. It lists the runtime per testbench and the simulation
status (``... ERROR``, ``FAILED``, ``NO ASSERTS`` or ``PASSED``). See
:doc:`Simulation </UsingPoC/Simulation>` for more details.


.. _QUICK:RunSynthesis:

Run a Synthesis
***************

The following quick example uses the Xilinx Systesis Tool (XST) to synthesize a
netlist for IP core ``arith_prng`` (Pseudo Random Number Generator - PRNG). The
VHDL file ``arith_prng.vhdl`` is located at ``PoCRoot\src\arith`` and virtually
a member in the `PoC.arith` namespace. So the module can be identified by an
unique name: ``PoC.arith.prng``, which is passed to the frontend script.

.. rubric:: Example:

.. code-block:: PowerShell

   cd PoCRoot
   .\poc.ps1 xst PoC.arith.prng --board=KC705

The CLI command ``xst`` chooses *Xilinx Synthesis Tool* as the synthesizer and
passes the fully qualified PoC entity name ``PoC.arith.prng`` as a parameter
to the tool. Additionally, the development board name is required to load the
correct ``my_config.vhdl`` file. All required source file are gathered and
synthesized to a netlist.

.. image:: /_static/images/xst/arith_prng.posh.png
   :target: /_static/images/xst/arith_prng.posh.png
   :alt: PowerShell console output after running PoC.arith.prng with XST.


.. _QUICK:Updating:

Updating
********

The PoC-Library can be updated by using ``git fetch`` and ``git merge``.

.. code-block:: PowerShell

   cd CaDHRoot
   # update the local repository
   git fetch --prune
   # review the commit tree and messages, using the 'treea' alias
   git treea
   # if all changes are OK, do a fast-forward merge
   git merge


.. seealso::
   :doc:`Running one or more testbenches </UsingPoC/Simulation>`
      The installation can be checked by running one or more of PoC's testbenches.
   :doc:`Running one or more netlist generation flows </UsingPoC/Synthesis>`
      The installation can also be checked by running one or more of PoC's
      synthesis flows.
