AIRFOIL=$1
POLARFILE=$1.ncrit09.t$2.Re$3k.pol

if [ -f $POLARFILE ] ; then
    rm $POLARFILE
fi

xfoil << EOF
load $1.dat
gdes
tset $2e-3


oper
v $3e3
pacc 
$POLARFILE

iter 200
aseq -3 11 1

quit
EOF

