RE="500 1000 1500 2000 2500 3000"
TAU="50 75 100 125 150 175"
FOIL=dai1336
for r in $RE
do
    for t in $TAU
    do
        ./genpolar.sh $FOIL $t $r
    done
done
