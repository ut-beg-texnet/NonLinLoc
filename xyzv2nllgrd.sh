#!/bin/bash
# xyz2nllgrd.sh
# Converts a velocity model expressed as an xyzv ascii file to the NLLoc grid format
# (which is composed by an ascii .hdr file and a binary .buf file)
# Requires a2b which is enclosed (not to be confused with the one from Seismic Unix)
#
# Set DOUBLE_PRECISION to 1 if you want double precision
#
# (c) 2010-2011 - Claudio Satriano <satriano@ipgp.fr>
export LC_ALL=C

DOUBLE_PRECISION=1

function getdelta() {
	awk '{
		if(old1) {
			delta=$1-old1
			if(olddelta) {
				ddelta=delta-olddelta;
				if(ddelta<0) ddelta=-ddelta;
				if(ddelta>0.01) {
					delta = "error"
					exit;	
				}
			}
			olddelta=delta
		}
		old1=$1
	} END {
		print delta
	}' $1
}

velocity=0
grid_type="SLOW_LEN"
if [ ""$2 == "-v" ]
then
	velocity=1
	grid_type="VELOCITY"
fi
if [ ""$1 == "-v" ]
then
	velocity=1
	grid_type="VELOCITY"
	shift
fi	
if [ -z ""$1 ]
then
	echo "Usage: $0 file.xyzv [-v]"
	echo " -v : output a velocity grid (default is slowness*lenght)"
	exit -1
fi
xyzvfile=$1
bname=`basename $xyzvfile .xyzv`

if [ "$DOUBLE_PRECISION" == "1" ]
then
	DP="-d"
fi

echo -n "Sorting file..."
sort $xyzvfile -nk1 -nk2 -nk3 > $bname.sort
echo " done"

# Compute delta, npts and origin for each direction
awk '{print $1}' $bname.sort | sort -nk1 | uniq > $bname.x 
nx=`cat $bname.x | wc | awk '{print $1}'`
origx=`cat $bname.x | head -n1 | awk '{print $1}'`
deltax=`getdelta $bname.x`
if [ ""$deltax == "error" ]
then
	echo "Error: deltax is not the same for all the samples"
	exit -1
fi
awk '{print $2}' $bname.sort | sort -nk1 | uniq > $bname.y 
ny=`cat $bname.y | wc | awk '{print $1}'`
origy=`cat $bname.y | head -n1 | awk '{print $1}'`
deltay=`getdelta $bname.y`
if [ ""$deltay == "error" ]
then
	echo "Error: deltay is not the same for all the samples"
	exit -1
fi
awk '{print $3}' $bname.sort | sort -nk1 | uniq > $bname.z 
nz=`cat $bname.z | wc | awk '{print $1}'`
origz=`cat $bname.z | head -n1 | awk '{print $1}'`
deltaz=`getdelta $bname.z`
if [ ""$deltaz == "error" ]
then
	echo "Error: deltaz is not the same for all the samples"
	exit -1
fi
# OK, done it! 

printf "%4d %4d %4d  %7.3f %7.3f %7.3f  %7.3f %7.3f %7.3f  %s\n" $nx $ny $nz $origx $origy $origz $deltax $deltay $deltaz  $grid_type > $bname.hdr

if [ $velocity -eq 1 ]
then
	awk '{print $4}' $bname.sort | ./a2b - $DP > $bname.buf
else
	# We use deltax to compute slowness*lenght. Is it correct???
	awk -v deltax=$deltax '{if($4) print (1/$4)*deltax; else print "NaN"}' $bname.sort | ./a2b - $DP > $bname.buf
fi

echo "Done! Output written to: $bname.hdr and $bname.buf"

rm -f $bname.sort $bname.[x,y,z]
