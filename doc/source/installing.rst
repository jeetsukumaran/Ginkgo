******************************************
Obtaining, Building, and Installing Ginkgo
******************************************

Obtaining Ginkgo
================

Ginkgo is available as a source distribution under the GNU General Public License (version 3) from the Git repository at:

    http://github.com/jeetsukumaran/Ginkgo

Source code archives of the latest releases can be found here:

    http://github.com/jeetsukumaran/Ginkgo/downloads

or the entire project repository can be cloned using Git via::

    $ git git://github.com/jeetsukumaran/Ginkgo.git

Building Ginkgo
===============

At the top level of the project directory, run::

    $ ./bootstrap.sh

and then::

    $ ./configure

By default, this will set up Ginkgo to be installed in its own sub-package directory in "``/opt/ginkgo``".
If you want to install the package in another location, you can specify this using the "``--prefix``" argument to configure.
For example, if user "alf" wishes to install the package in the subdirectory "``programs``" of his home directory::

    $ ./configure --prefix=/home/alf/programs/

After that, build the package::

    $ make

Installing Ginkgo
=================

To install the Ginkgo package, run::

    $ make install

If the directory in which Ginkgo is to be installed is not owned by the user, then the following is neccessary instead::

    $ sudo make install

Note that the ``ginkgo`` executable, which is the actual simulator itself, will most likely *not* be on the system path, as it is located in the "``ginkgo/bin``" subdirectory of the install location.
For example, if the installation prefix was "``/opt/ginkgo``", then the full path of the executable will be "``/opt/ginkgo/bin/ginkgo``".
So, for user "alf" above, the path to the executable will be "``/home/alf/programs/ginkgo/bin/ginkgo``".
For convenience, you can modify your ``$PATH`` environmental variable to include the ginkgo installation binary directory, or simply copy the executable to your system path, e.g.::

    $ sudo cp /opt/ginkgo/bin/ginkgo /usr/local/bin


