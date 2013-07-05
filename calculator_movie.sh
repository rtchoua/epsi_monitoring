#!/bin/bash
#
# Create FLV movies from a set of png/jpg images using mencoder, ffmpeg 
#
# Call with the directory that contains the images. All variables under this
# directory will be processed and AVI and FLV generated from PNG or JPG images
#
# Author: Norbert Podhorszki, pnorbert@ornl.gov, 2008
#


umask 022

######################
# set default values #
######################
KEEPAVI=yes
FPS=16
KEYINT=1
#LOG=/dev/stdout
IMGDIR=.
IMGTYPE=png
_SCRIPTDIR_=`dirname $0`            # relative (or absolute) path of this script
SCRIPTDIR=`cd $_SCRIPTDIR_; pwd`    # absolute path of this script



#############
# Functions #
#############

#
# Usage printout
#
function Usage () {
    echo "Usage: `basename $0` [-q] [-a] [-r framerate] [-k keyint] [-d dir] [-t imgtype] [-l log] varname"
    echo "Arguments"
    echo "  -d dir        The directory which contains the images (default: $IMGDIR)."
    echo "  -r framerate  Frame/sec speed for the movie. Default: $FPS fps."
    echo "  -k keyint     Interval between keyframes. Default is $KEYINT."
    echo "                   Try to choose a divisor of n-1 whenever possible,"
    echo "                   where n = the number of images to be encoded to FLV, "
    echo "                   otherwise the last frame of the FLV movie will not be seekable."
    echo "  -a            Remove the intermediate AVI movie. It is kept by default."
    echo "  -t imgtype    Type of img (png or jpg). Default is $IMGTYPE."
    echo "  -l log        Log file, default: stdout. "
    echo "                   Note that mencode, ffmpeg and flvxmltag makes separate log"
    echo "                   files in the -d dir."
    echo "  -q            Quiet mode (no logs)"
    echo "  -h            This help"
}

#
# add a path to $PATH
#
function pathmunge () {
        if ! echo $PATH | /bin/egrep -q "(^|:)$1($|:)" ; then
           if [ "$2" = "after" ] ; then
              PATH=$PATH:$1
           else
              PATH=$1:$PATH
           fi
        fi
}

#
# add a path to LD_LIBRARY_PATH
#
function ldmunge () {
        if ! echo $LD_LIBRARY_PATH | /bin/egrep -q "(^|:)$1($|:)" ; then
           if [ "$2" = "after" ] ; then
              LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$1
           else
              LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
           fi
        fi
}

#####################
# Process arguments #
#####################

# process option arguments
while getopts ":d:r:k:a:l:i:t:qh" Option
do
  case $Option in
        d) IMGDIR=$OPTARG;;
        r) FPS=$OPTARG;;
        k) KEYINT=$OPTARG;;
        a) KEEPAVI=no;;
        l) LOG=$OPTARG;;
        t) IMGTYPE=$OPTARG;;
        q) VERBOSE=no;LOG=/dev/null;;
        h) Usage; exit 0;;
        *) echo "Invalid option -$Option."; Usage; exit 255;;   # DEFAULT
  esac
done
shift $(($OPTIND - 1))

[ -z "$1" ] && echo "Give the variable name" 1>&2 && Usage && exit 254
VARNAME=$1
LOG=$VARNAME.log

############################
# Prepare path and ld path #
############################

. /etc/profile.d/modules.sh
#module load flvtool2
#module load ffmpeg
#module load mplayer

pathmunge /sw/sith/flvtool2/1.0.6/centos5.5_gnu4.1.2/bin/flvtool2/flvtool2/bin
pathmunge /sw/sith/mplayer/1.0rc3/centos5.5_gnu4.4.4/bin 
ldmunge /sw/sith/flvtool2/1.0.6/centos5.5_gnu4.1.2/bin/flvtool2/flvtool2/lib
ldmunge /sw/sith/mplayer/1.0rc3/centos5.5_gnu4.4.4/bin
#ldmunge /apps/ffmpeg/lame/lib
#ldmunge /apps/ffmpeg/lib
#ldmunge /data/web2/dev/rbarreto/lib


################
# Create movie #
################

MYEX=0
    
# enter the directory
pushd $IMGDIR >/dev/null

echo "This is NORBERT's script $0" > $LOG
echo "Path is $PATH" >> $LOG
echo "LDPath is $LD_LIBRARY_PATH" >> $LOG
#echo "ldd ffmpeg:" >> $LOG
#ldd `which ffmpeg` >> $LOG

# create the AVI movie from the images
AVI=$VARNAME.avi
AVI_KEYINT=160
AVI_VBITRATE=512000
echo "/sw/sith/mplayer/1.0rc3/centos5.5_gnu4.4.4/bin/mencoder "mf://$VARNAME*.$IMGTYPE" -mf type=$IMGTYPE:fps=$FPS -ovc lavc -lavcopts vcodec=mpeg4:autoaspect=1:vbitrate=$AVI_VBITRATE:mbd=2:keyint=$AVI_KEYINT -nosound -o $AVI" >> $LOG
/sw/sith/mplayer/1.0rc3/centos5.5_gnu4.4.4/bin/mencoder "mf://$VARNAME*.$IMGTYPE" -mf type=$IMGTYPE:fps=$FPS -ovc lavc -lavcopts vcodec=mpeg4:autoaspect=1:vbitrate=$AVI_VBITRATE:mbd=2:keyint=$AVI_KEYINT -nosound -o $AVI &>$VARNAME.log.mencoder
EX=$?
echo "Exit code of mencoder = $EX" >> $LOG
[ $EX != 0 ] && MYEX=$EX

# create FLV from AVI
FLV=$VARNAME.flv
FLV_BITRATE=512k
echo "/sw/sith/ffmpeg/0.6.1/centos5.5_gnu4.4.4/bin/ffmpeg -i $AVI -b $FLV_BITRATE -g $KEYINT -r $FPS -y $FLV" >> $LOG
/ccs/proj/e2e/rbarreto/ffmpeg/ffmpeg -i $AVI -b $FLV_BITRATE -g $KEYINT -r $FPS -y $FLV &>$VARNAME.log.ffmpeg
EX=$?
echo "Exit code of ffmpeg = $EX" >> $LOG
[ $EX != 0 ] && MYEX=$EX

# remove the intermediate AVI if not asked to keep it
[ "$KEEPAVI" == "no" ] && rm -f $AVI

# Tag the flv movie
TAGSCRIPT=$SCRIPTDIR/flvxmltags.py
if [ -f $TAGSCRIPT ]; then
    echo "python $TAGSCRIPT --input=$FLV" >>$LOG
    EX=$?
    python $TAGSCRIPT --input=$FLV  &>$VARNAME.log.flvxmltags
    echo "Exit code of flvxmltags = $EX" >> $LOG
    [ $EX != 0 ] && MYEX=$EX
else
    echo "Script $TAGSCRIPT not found." >>$LOG
fi

# return to original dir 
popd >/dev/null

exit $MYEX
