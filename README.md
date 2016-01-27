#checkcel

The checkcel program allows rapid validity testing of Affymetrix `.CEL` files. The aim of this program is to bulk validate multiple `.CEL` files from public data.

##Licence

checkcel is released under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl.html).

##Options & Arguments

checkcel is called as follows:

    checkcel [-cfvh] file [...]

* `-h`: print help
* `-v`: print version
* `-f`: filter out invalid `.CEL` files
* `-c`: calculate & display intensity statistics

##Output Format

Each output line gives data for a single input file. Output data are tab-delimited. If `-f` is specified, only valid files are returned, otherwise invalid files are returned with the value `unknown`. If `-c` is not specified, the output columns are:

* `CEL` file name
* file format
* chip ID
* creation algorithm
* row count
* column count
* margin
* outier cell count
* masked cell count

If `-c` is specified, the following columns are appended to each line:

* minimum intensity value
* maximum intensity value
* unique value count
* invalid value count

##Building checkcel

checkcel should be made by:

    cd checkcel
    ./make

