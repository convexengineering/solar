AIRFOIL=$1
POLARFILE=$1.ncrit09.Re$2k.pol

if [ -f $POLARFILE ] ; then
    rm $POLARFILE
fi

xfoil << EOF
load $1.dat
oper
v $2e3
pacc 
$POLARFILE

iter 200
aseq 5 11 0.5

quit
EOF

