******************************************
Obtaining, Building, and Installing Ginkgo
******************************************

Obtaining and Installing Pre-Built Ginkgo Binaries
==================================================

Pre-built Ginkgo binaries (executables) are available for the Apple OS X and Microsoft Windows platforms:

    * :download:`Apple OS X (32-bit) </_downloads/ginkgo-3.9.0-bin-osx32.zip>`
    * :download:`Apple OS X (64-bit) </_downloads/ginkgo-3.9.0-bin-osx64.zip>`
    * :download:`Microsoft Windows </_downloads/ginkgo-3.9.0-bin-w32.zip>`

If you are running "Snow Leopard" (OS X 10.6) or higher, you want the :download:`64-bit OS X version </_downloads/ginkgo-3.9.0-bin-osx64.zip>`, otherwise you want the :download:`32-bit version </_downloads/ginkgo-3.9.0-bin-osx32.zip>`.

As is typical for most other platforms, such as the various flavors of Linux, you would need to build and install directly from the source code.
See below for details on how to :ref:`obtain, build and install Ginkgo directly from the source code <ginkgo-from-source>`.

"Installing" the pre-built binaries involves nothing more than downloading the appropriate archive, unzipping the program, and copying the program to a suitable location.
On the OS X platforms, a "suitable location" is typically "``/usr/local/bin``", so installing Ginkgo would entail something like the following::

    $ curl -O http://www.phylo.bio.ku.edu/ginkgo/_downloads/ginkgo-3.9.0-bin-osx64.zip
    $ unzip ginkgo-3.9.0-bin-osx64.zip
    $ sudo mv ginkgo /usr/local/bin

The "``sudo``" is necessary because you are copying the file to a system-wide directory.
You do not actually need to do this, as Ginkgo runs quite well from any location, but by doing this you avoid the need to enter the full path (i.e., parent directories) when invoking Ginkgo.

.. _ginkgo-from-source:

Obtaining, Building and Installing Ginkgo from the Source Code
==============================================================

Downloading the Source Code
---------------------------

The Ginkgo source code is distributed under the GNU General Public License (version 3).

You can download the snapshot archive of the latest public release of the source code from the Git repository at:

    http://github.com/jeetsukumaran/Ginkgo

Source code archives of the latest releases can be found here:

    http://github.com/jeetsukumaran/Ginkgo/downloads

or the entire project repository can be cloned using <a href="http://www.git-scm.com">Git<a> via the command::

    $ git git://github.com/jeetsukumaran/Ginkgo.git

Building Ginkgo
---------------

At the top level of the project directory, run::

    $ ./bootstrap.sh

and then::

    $ ./configure

By default, this will set up Ginkgo to be installed in the "``bin``" sub-directory of "``/usr/local``".
If you want to install the program in another location, you can specify this using the "``--prefix``" argument to configure.
For example, if user "alf" wishes to install the package in the "``bin``" subdirectory of "``programs``" in his home directory::

    $ ./configure --prefix=/home/alf/programs/

After that, build the package::

    $ make

Installing Ginkgo
-----------------

To install the Ginkgo package, run::

    $ make install

If the directory in which Ginkgo is to be installed is not owned by the user, then the following necessary instead::

    $ sudo make install

As noted above, by default (if the "``--prefix``" argument was not passed to "``configure``" above), this will result in the ``ginkgo`` executable being placed in the "``/usr/local/bin``" directory.
In most typical system configurations, this will allow for the program to be available system-wide for all users, but it will also require that you have administrative writing privileges (hence the "``sudo``" in the command above).

If you do not want to install ``ginkgo`` on a system-wide path, or you do not not have administrative privileges, then you will have to pass the installation prefix to the "``configure``" command as described above.
In these cases, you probably want to modify your ``$PATH`` environmental variable to include the ginkgo installation binary directory.

Windows Users
-------------

None of the above really applies to Windows users.
If you are feeling adventurous, you could trying creating a Visual Studio project, importing all the files in the "``ginkgocc``" source directory, and building the result.
It should work --- that is how the Windows binaries linked to on this page were built.
However, it probably makes more sense to simply use one of the pre-built binaries.
