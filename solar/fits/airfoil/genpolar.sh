AIRFOIL=$1
REN=$2
NCRIT=$3
POLARFILE=$1.ncrit0$3.Re$2k.pol

if [ -f $POLARFILE ] ; then
    rm $POLARFILE
fi

xfoil << EOF
load $1.dat
oper
v $2e3
vpar
n $3

pacc 
$POLARFILE

iter 200
aseq 5 11 0.5

quit
EOF

