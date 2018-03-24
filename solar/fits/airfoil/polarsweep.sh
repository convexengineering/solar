RE="125 150 200 300 400 500 600"
TAU="80 90 100 105"
FOIL=dai1336a
for r in $RE
do
    for t in $TAU
    do
        ./genpolar.sh $FOIL.$t $r
    done
done
