#!/bin/tcsh
setenv TZ UTC
set BIN = $HOME/src
cd $BIN

setenv PGMSCALER 1000
# setenv GENERATE_INTER_ACC True

set TIMESTAMP = $1
# setenv TIMESTEPS 18
set INTSTEPS = 10
set PPNDIR = /dev/shm/ppn
set OUTDIR = $PPNDIR/acc
set PROD = fmi.radar.composite.lowest_FIN_RAVAKE.pgm
set OBSDIR = /mnt/meru/data/prod/radman/latest/fmi/radar/composite/lowest
set OBSFILE = $OBSDIR/"$TIMESTAMP"_"$PROD"
set NCFILE = $PPNDIR/nc_"$TIMESTAMP".h5
set MEMBERS = `h5dump -a /meta/configuration/ENSEMBLE_SIZE $NCFILE | grep : | cut -d: -f2 | tr -d '"'`
set LASTMEMBER = `expr $MEMBERS \- 1`
set PIDSEEK = "$TIMESTAMP"_'member-*'.pid

if(! -e $OUTDIR) then 
   mkdir $OUTDIR
else
   find $OUTDIR/ -name '*.???' -mmin +15 -exec rm -f {} \;
endif

echo '========='
date -u

echo "Starting interpolation processes for $TIMESTAMP with $MEMBERS members\n" 

set M = 0
while($M <= $LASTMEMBER)
   set MM = $M
   if($M < 10) set MM = 0"$M"
   set MEMBERNAME = member-"$MM"
   echo "Starting $MEMBERNAME"
   set EXEC = "time ./member_interp $TIMESTAMP $MEMBERNAME $NCFILE $OBSFILE $OUTDIR $INTSTEPS"
   if($M < $LASTMEMBER) then
       # background processes for members except the last one
       $EXEC > "$MEMBERNAME".log & 
   endif
   set M = `expr $M + 1`
end

# foreground process for the last member
$EXEC > "$MEMBERNAME".log 

echo "\nLast member finished, waiting others (if any) to finish:"

set FAIL = 0
set WAITS = 60
while(1)
   set PIDFILES = `ls $PIDSEEK`
   if("$PIDFILES" == "") break
   set PIDS = `cat $PIDFILES`
   echo "Still running:" 
   echo "$PIDFILES" | tr ' ' "\n"
   sleep 0.5
   set WAITS = `expr $WAITS \- 1`
   if(! $WAITS) then
      echo "Killing jammed processes $PIDS"
      kill -KILL $PIDS
      echo "Interpolation $TIMESTAMP failed!" >> interp_err.log
      exit 1
   endif
end

echo "\nInterpolations finished."

echo "\nConcatenation of member accumulation files:"

cd $OUTDIR
foreach LF (member-00_RAVACC_"$TIMESTAMP"*.dat)
   set ACCFILE = `echo $LF | cut -d_ -f2-`
   cat member-??_"$ACCFILE" > $ACCFILE
   echo "$ACCFILE ready"
end
rm member-*.dat

echo "\nInterpolation procedure of $TIMESTAMP done"
date -u

exit

