--------------
About HiC Tools
--------------
HiC Tools is a software package used to build `.hic` files. This distribution includes the source code for this set of
tools and builds the latest `v9` `.hic` files.

Previously, this software used to be bundled with Juicer Tools in the Juicebox repo.

HiC Tools was originally created by <a href="https://github.com/jrobinso">Jim Robinson</a>,
<a href="https://github.com/nchernia">Neva C. Durand</a>, and <a href="http://www.erez.com/">Erez Lieberman Aiden</a>.
It has been substantially upgraded and improved by <a href="https://github.com/sa501428">Muhammad S Shamim</a>,
<a href="https://github.com/suhas-rao">Suhas S.P. Rao</a>, and many additional contributors.

**If you use HiC Tools in your research, please cite:
Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman
Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments."
Cell Systems 3(1), 2016.**


--------------
Questions?
--------------

For FAQs, or for asking new questions, please see our forum: <a href="https://aidenlab.org/forum.html">aidenlab.org/forum.html</a>.

--------------
IntelliJ Setup
--------------

Use IntelliJ IDEA (Community edition - free)

To set up in IDEA, have the Java SDK installed
then you'll point to it (IntelliJ has lots of documentation on this sort of thing).

* Then go to `VCS` -> `checkout from version control`.
* Then go to `Run` -> `Edit Configurations`.
* With the `+` sign, add `Application`.
* You'll create one of these for HiC Tools.
* Set the main class by clicking the little `...` button next to the text box for main class
* HiCTools.java is the main method class.
* Under VM Options, set `-Xmx2000m`
* The `Xmx2000m` flag sets the maximum memory heap size to 2GB. Depending on your computer you might want more or less.
  Some tools will break if there's not enough memory and the file is too large, but don't worry about that for
  development; 2GB should be fine.
* One last note: be sure to `Commit and Push` when you commit files, it's hidden in the dropdown menu button in the
  commit window.

----------------------------------
Hardware and Software Requirements
----------------------------------
The minimum software requirement to run Juicebox is a working Java installation
(version > 1.8) on Windows, Linux, and Mac OSX. We recommend using the latest Java version available, but please do not
use the Java Beta Version. Minimum system requirements for running Java can be found at
https://java.com/en/download/help/sysreq.xml. To download and install the latest Java Runtime Environment (JRE), please
go to https://www.java.com/download.

-------------
Documentation
-------------
We have extensive documentation for how to use HiC Tools at
https://github.com/aidenlab/juicer/wiki.
